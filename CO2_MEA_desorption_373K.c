/*============================================================
 * File   : co2_mea_desorption_udf.c
 * Brief  : UDF for CO2 desorption from MEA solution
 *          in a gas-liquid multiphase flow (Euler-Euler)
 *
 * Species index:
 *   Liquid phase  -- MEA (0), MEAH (1), MEACOO- (2), H2O (2)
 *   Gas    phase  -- CO2 (0), N2  (1)
 *
 * UDMI storage:
 *   UDMI[0] -- liquid dynamic viscosity      [Pa·s]
 *   UDMI[1] -- turbulent energy dissipation  [m²/s³]
 *   UDMI[2] -- CO2 mass source term          [kg/(m³·s)]
 *============================================================*/

#include "udf.h"
#include "mem.h"
#include "sg_mphase.h"
#include "metric.h"
#include "global.h"
#define density 1050    /* liquid phase density                         [kg/m³] */
#define K       0.167   /* coefficient in small eddy model      */

/*------------------------------------------------------------
 * DEFINE_PROPERTY : liq_mu_udf
 * Returns : liquid dynamic viscosity  [Pa·s]
 * Method  : empirical correlation as a function of
 *           liquid temperature T_l and CO2 loading
 *------------------------------------------------------------*/

DEFINE_PROPERTY(liq_mu_udf, cell, thread_liq)
{
    real T_l, loading, v;
 
    T_l           = C_T(cell, thread_liq);               /* liquid temperature           [K]  */
    real y_meacoo = C_YI(cell, thread_liq, 2);           /* MEACOO- mass fraction        [-]  */
    real denom    = 104.08 - 44.01 * y_meacoo;           /* denominator for molar mass        */ 
    loading = y_meacoo * 203.6 / denom;                  /* CO2 loading  [mol CO2 / mol MEA]  */
    v       = exp(7.08 - 0.021 * T_l + 0.97 * loading) / 1000.0;   /* dynamic viscosity [Pa·s] */
 
    C_UDMI(cell, thread_liq, 0) = v;                     /* store to UDMI[0] for post-processing */
    return v;
}

/*------------------------------------------------------------
 * DEFINE_SOURCE : CO2_src
 * Returns : CO2 mass source term in the GAS phase  [kg/(m³·s)]
 *           Positive -> gas gains CO2 (desorption)
 *           Negative -> gas loses CO2 (absorption)
 *------------------------------------------------------------*/
 DEFINE_SOURCE(CO2_src, cell, thread_gas, dS, eqn)
{
    Thread *thread_mix, *thread_liq;
    thread_mix = THREAD_SUPER_THREAD(thread_gas);        /* mixture thread              */
    thread_liq = THREAD_SUB_THREAD(thread_mix, 0);       /* liquid sub-thread (index 0) */
 
    real ae, c_mea, d_n2o, d_co2, c_d, kL, k_r, Ha, E;
    real yL, m_dot_co2, Henry, loading, v, c_meacoo;
 
    real T_l = C_T(cell, thread_liq);                    /* liquid temperature  [K] */
	if (C_VOF(cell, thread_gas) < 0.99 && C_VOF(cell, thread_gas) > 0.01)
    {
        /* --- interfacial area density (spherical bubble assumption) --- */
        /* ae = 6 * alpha_g / d_b  [m²/m³] */
        ae = 6.0 * C_VOF(cell, thread_gas) / C_PHASE_DIAMETER(cell, thread_gas);
		
		 /* --- species and CO2 loading --- */
        real y_meacoo    = C_YI(cell, thread_liq, 2);          /* MEACOO- mass fraction  [-]     */
        real denom = 104.08 - 44.01 * y_meacoo;
        loading    = y_meacoo * 203.6 / denom;                 /* CO2 loading  [mol/mol]         */
		
		/* --- liquid dynamic viscosity --- */
        v = exp(7.08 - 0.021 * T_l + 0.97 * loading) / 1000.0;     /* [Pa·s] */
		
		/* --- molar concentrations --- */
        c_mea    = C_YI(cell, thread_liq, 0) * density / 61.0;      /* MEA     [mol/L] */
        c_meacoo = C_YI(cell, thread_liq, 2) * density / 104.0;     /* MEACOO- [mol/L] */
		
		/* --- diffusivity via N2O analogy --- */
        d_n2o = (5.07 + 0.865 * c_mea + 0.278 * c_mea * c_mea)
                * exp((-2371.0 - 93.4 * c_mea) / T_l) * 1e-6;      /* D_N2O  [m²/s] */
        d_co2 = d_n2o * 0.4635 * exp(252.0 / T_l);                  /* D_CO2  [m²/s] */
		
		/* --- turbulent energy dissipation rate --- */
        /* epsilon = C_mu * k * omega  [m²/s³] */
        c_d = 0.09 * C_K(cell, thread_liq) * C_O(cell, thread_liq);
		
		/* --- liquid-side mass transfer coefficient (Calderbank-Moo-Young) --- */
        kL = K * pow(d_co2, 0.5) * pow((c_d / v * density), 0.25); /* [m/s] */
		
		 /* --- pseudo-first-order reaction rate constant --- */
        /* k_r = k2 * c_MEA  [1/s] */
        k_r = pow(10.0, (10.99 - 2152.0 / T_l)) * c_mea;
		
		/* --- Hatta number and enhancement factor --- */
        Ha = pow((k_r * d_co2), 0.5) / kL;              /* Hatta number          [-] */
        E  = pow((1.0 + pow(Ha, 2)), 0.5);               /* enhancement factor (fast reaction approximation) */
		
		/* --- Henry's constant --- */
        Henry = exp(896.5 / T_l - 10.5913);              /* [mol/(m³·kPa)] */
		
		/* --- driving force yL = equilibrium CO2 conc. - actual CO2 conc.  [mol/m³] --- */
        yL = exp(35.81 * y_meacoo - 3.09) * Henry
             - (C_P(cell, thread_gas) / 1000.0 + 101.325) * Henry
             * (C_YI(cell, thread_gas, 0)
             / (C_YI(cell, thread_gas, 0) + 2.44 * (1.0 - C_YI(cell, thread_gas, 0))));
			 
		 /* --- CO2 gas-phase mass source term --- */
        /* E * kL [m/s] * ae [m²/m³] * M_CO2 [kg/mol] * yL [mol/m³] = [kg/(m³·s)] */
        m_dot_co2 = E * kL * ae * 44.0 * yL;
		
		   dS[eqn] = 0.0;                
    }
    else
    {
        /* source set to zero when VOF is outside the valid interface range */
		c_d = 0.0;
        m_dot_co2 = 0.0;
        dS[eqn]   = 0.0;
    }
	C_UDMI(cell, thread_liq, 1) = c_d;
    C_UDMI(cell, thread_liq, 2) = m_dot_co2; 
    return m_dot_co2;
}
		
	
/*------------------------------------------------------------
 * DEFINE_SOURCE : gas_src
 * Returns : CO2 mass source term in the GAS phase  [kg/(m³·s)]
 *           Positive -> gas gains CO2 (desorption)
 *           Negative -> gas loses CO2 (absorption)
 *------------------------------------------------------------*/
 DEFINE_SOURCE(gas_src, cell, thread_gas, dS, eqn)
{
    Thread *thread_mix, *thread_liq;
    thread_mix = THREAD_SUPER_THREAD(thread_gas);        /* mixture thread              */
    thread_liq = THREAD_SUB_THREAD(thread_mix, 0);       /* liquid sub-thread (index 0) */
 
    real ae, c_mea, d_n2o, d_co2, c_d, kL, k_r, Ha, E;
    real yL, m_dot_co2, Henry, loading, v, c_meacoo;
 
    real T_l = C_T(cell, thread_liq);                    /* liquid temperature  [K] */
	if (C_VOF(cell, thread_gas) < 0.99 && C_VOF(cell, thread_gas) > 0.01)
    {
        /* --- interfacial area density (spherical bubble assumption) --- */
        /* ae = 6 * alpha_g / d_b  [m²/m³] */
        ae = 6.0 * C_VOF(cell, thread_gas) / C_PHASE_DIAMETER(cell, thread_gas);
		
		 /* --- species and CO2 loading --- */
        real y_meacoo    = C_YI(cell, thread_liq, 2);          /* MEACOO- mass fraction  [-]     */
        real denom = 104.08 - 44.01 * y_meacoo;
        loading    = y_meacoo * 203.6 / denom;                 /* CO2 loading  [mol/mol]         */
		
		/* --- liquid dynamic viscosity --- */
        v = exp(7.08 - 0.021 * T_l + 0.97 * loading) / 1000.0;     /* [Pa·s] */
		
		/* --- molar concentrations --- */
        c_mea    = C_YI(cell, thread_liq, 0) * density / 61.0;      /* MEA     [mol/L] */
        c_meacoo = C_YI(cell, thread_liq, 2) * density / 104.0;     /* MEACOO- [mol/L] */
		
		/* --- diffusivity via N2O analogy --- */
        d_n2o = (5.07 + 0.865 * c_mea + 0.278 * c_mea * c_mea)
                * exp((-2371.0 - 93.4 * c_mea) / T_l) * 1e-6;      /* D_N2O  [m²/s] */
        d_co2 = d_n2o * 0.4635 * exp(252.0 / T_l);                  /* D_CO2  [m²/s] */
		
		/* --- turbulent energy dissipation rate --- */
        /* epsilon = C_mu * k * omega  [m²/s³] */
        c_d = 0.09 * C_K(cell, thread_liq) * C_O(cell, thread_liq);
		
		/* --- liquid-side mass transfer coefficient (Calderbank-Moo-Young) --- */
        kL = K * pow(d_co2, 0.5) * pow((c_d / v * density), 0.25); /* [m/s] */
		
		 /* --- pseudo-first-order reaction rate constant --- */
        /* k_r = k2 * c_MEA  [1/s] */
        k_r = pow(10.0, (10.99 - 2152.0 / T_l)) * c_mea;
		
		/* --- Hatta number and enhancement factor --- */
        Ha = pow((k_r * d_co2), 0.5) / kL;              /* Hatta number          [-] */
        E  = pow((1.0 + pow(Ha, 2)), 0.5);               /* enhancement factor (fast reaction approximation) */
		
		/* --- Henry's constant --- */
        Henry = exp(896.5 / T_l - 10.5913);              /* [mol/(m³·kPa)] */
		
		/* --- driving force yL = equilibrium CO2 conc. - actual CO2 conc.  [mol/m³] --- */
        yL = exp(35.81 * y_meacoo - 3.09) * Henry
             - (C_P(cell, thread_gas) / 1000.0 + 101.325) * Henry
             * (C_YI(cell, thread_gas, 0)
             / (C_YI(cell, thread_gas, 0) + 2.44 * (1.0 - C_YI(cell, thread_gas, 0))));
			 
		 /* --- CO2 gas-phase mass source term --- */
        /* E * kL [m/s] * ae [m²/m³] * M_CO2 [kg/mol] * yL [mol/m³] = [kg/(m³·s)] */
        m_dot_co2 = E * kL * ae * 44.0 * yL;
		
		   dS[eqn] = 0.0;                
    }
    else
    {
        /* source set to zero when VOF is outside the valid interface range */
        m_dot_co2 = 0.0;
        dS[eqn]   = 0.0;
    }

    return m_dot_co2;
}	

/*------------------------------------------------------------
 * DEFINE_SOURCE : mea_src
 * Returns : mea mass source term in the Liq phase  [kg/(m³·s)]
 *           Positive -> liquid generates MEA (desorption)
 *           Negative -> liquid consumes MEA (absorption)
 *------------------------------------------------------------*/
 DEFINE_SOURCE(liqmea_src, cell, thread_liq, dS, eqn)
{
    Thread *thread_mix, *thread_gas;
    thread_mix = THREAD_SUPER_THREAD(thread_liq);        /* mixture thread              */
    thread_gas = THREAD_SUB_THREAD(thread_mix, 1);       /* gas sub-thread (index 1) */
 
    real ae, c_mea, d_n2o, d_co2, c_d, kL, k_r, Ha, E;
    real yL, m_dot_mea, Henry, loading, v, c_meacoo;
 
    real T_l = C_T(cell, thread_liq);                    /* liquid temperature  [K] */
	if (C_VOF(cell, thread_gas) < 0.99 && C_VOF(cell, thread_gas) > 0.01)
    {
        /* --- interfacial area density (spherical bubble assumption) --- */
        /* ae = 6 * alpha_g / d_b  [m²/m³] */
        ae = 6.0 * C_VOF(cell, thread_gas) / C_PHASE_DIAMETER(cell, thread_gas);
		
		 /* --- species and CO2 loading --- */
        real y_meacoo    = C_YI(cell, thread_liq, 2);          /* MEACOO- mass fraction  [-]     */
        real denom = 104.08 - 44.01 * y_meacoo;
        loading    = y_meacoo * 203.6 / denom;                 /* CO2 loading  [mol/mol]         */
		
		/* --- liquid dynamic viscosity --- */
        v = exp(7.08 - 0.021 * T_l + 0.97 * loading) / 1000.0;     /* [Pa·s] */
		
		/* --- molar concentrations --- */
        c_mea    = C_YI(cell, thread_liq, 0) * density / 61.0;      /* MEA     [mol/L] */
        c_meacoo = C_YI(cell, thread_liq, 2) * density / 104.0;     /* MEACOO- [mol/L] */
		
		/* --- diffusivity via N2O analogy --- */
        d_n2o = (5.07 + 0.865 * c_mea + 0.278 * c_mea * c_mea)
                * exp((-2371.0 - 93.4 * c_mea) / T_l) * 1e-6;      /* D_N2O  [m²/s] */
        d_co2 = d_n2o * 0.4635 * exp(252.0 / T_l);                  /* D_CO2  [m²/s] */
		
		/* --- turbulent energy dissipation rate --- */
        /* epsilon = C_mu * k * omega  [m²/s³] */
        c_d = 0.09 * C_K(cell, thread_liq) * C_O(cell, thread_liq);
		
		/* --- liquid-side mass transfer coefficient (Calderbank-Moo-Young) --- */
        kL = K * pow(d_co2, 0.5) * pow((c_d / v * density), 0.25); /* [m/s] */
		
		 /* --- pseudo-first-order reaction rate constant --- */
        /* k_r = k2 * c_MEA  [1/s] */
        k_r = pow(10.0, (10.99 - 2152.0 / T_l)) * c_mea;
		
		/* --- Hatta number and enhancement factor --- */
        Ha = pow((k_r * d_co2), 0.5) / kL;              /* Hatta number          [-] */
        E  = pow((1.0 + pow(Ha, 2)), 0.5);               /* enhancement factor (fast reaction approximation) */
		
		/* --- Henry's constant --- */
        Henry = exp(896.5 / T_l - 10.5913);              /* [mol/(m³·kPa)] */
		
		/* --- driving force yL = equilibrium CO2 conc. - actual CO2 conc.  [mol/m³] --- */
        yL = exp(35.81 * y_meacoo - 3.09) * Henry
             - (C_P(cell, thread_gas) / 1000.0 + 101.325) * Henry
             * (C_YI(cell, thread_gas, 0)
             / (C_YI(cell, thread_gas, 0) + 2.44 * (1.0 - C_YI(cell, thread_gas, 0))));
			 
		 /* --- MEA liq-phase mass source term --- */
        /* 2 * E * kL [m/s] * ae [m²/m³] * M_MEA [kg/mol] * yL [mol/m³] = [kg/(m³·s)] */
        m_dot_mea = 2* E * kL * ae * 61.0 * yL;
		
		   dS[eqn] = 0.0;                
    }
    else
    {
        /* source set to zero when VOF is outside the valid interface range */
        m_dot_mea = 0.0;
        dS[eqn]   = 0.0;
    }

    return m_dot_mea;
}	

/*------------------------------------------------------------
 * DEFINE_SOURCE : meah_src
 * Returns : meah mass source term in the Liq phase  [kg/(m³·s)]
 *           Positive -> liquid generates MEAH (absorption)
 *           Negative -> liquid consumes MEAH (desorption)
 *------------------------------------------------------------*/
 DEFINE_SOURCE(liqmeah_src, cell, thread_liq, dS, eqn)
{
    Thread *thread_mix, *thread_gas;
    thread_mix = THREAD_SUPER_THREAD(thread_liq);        /* mixture thread              */
    thread_gas = THREAD_SUB_THREAD(thread_mix, 1);       /* gas sub-thread (index 1) */
 
    real ae, c_mea, d_n2o, d_co2, c_d, kL, k_r, Ha, E;
    real yL, m_dot_meah, Henry, loading, v, c_meacoo;
 
    real T_l = C_T(cell, thread_liq);                    /* liquid temperature  [K] */
	if (C_VOF(cell, thread_gas) < 0.99 && C_VOF(cell, thread_gas) > 0.01)
    {
        /* --- interfacial area density (spherical bubble assumption) --- */
        /* ae = 6 * alpha_g / d_b  [m²/m³] */
        ae = 6.0 * C_VOF(cell, thread_gas) / C_PHASE_DIAMETER(cell, thread_gas);
		
		 /* --- species and CO2 loading --- */
        real y_meacoo    = C_YI(cell, thread_liq, 2);          /* MEACOO- mass fraction  [-]     */
        real denom = 104.08 - 44.01 * y_meacoo;
        loading    = y_meacoo * 203.6 / denom;                 /* CO2 loading  [mol/mol]         */
		
		/* --- liquid dynamic viscosity --- */
        v = exp(7.08 - 0.021 * T_l + 0.97 * loading) / 1000.0;     /* [Pa·s] */
		
		/* --- molar concentrations --- */
        c_mea    = C_YI(cell, thread_liq, 0) * density / 61.0;      /* MEA     [mol/L] */
        c_meacoo = C_YI(cell, thread_liq, 2) * density / 104.0;     /* MEACOO- [mol/L] */
		
		/* --- diffusivity via N2O analogy --- */
        d_n2o = (5.07 + 0.865 * c_mea + 0.278 * c_mea * c_mea)
                * exp((-2371.0 - 93.4 * c_mea) / T_l) * 1e-6;      /* D_N2O  [m²/s] */
        d_co2 = d_n2o * 0.4635 * exp(252.0 / T_l);                  /* D_CO2  [m²/s] */
		
		/* --- turbulent energy dissipation rate --- */
        /* epsilon = C_mu * k * omega  [m²/s³] */
        c_d = 0.09 * C_K(cell, thread_liq) * C_O(cell, thread_liq);
		
		/* --- liquid-side mass transfer coefficient (Calderbank-Moo-Young) --- */
        kL = K * pow(d_co2, 0.5) * pow((c_d / v * density), 0.25); /* [m/s] */
		
		 /* --- pseudo-first-order reaction rate constant --- */
        /* k_r = k2 * c_MEA  [1/s] */
        k_r = pow(10.0, (10.99 - 2152.0 / T_l)) * c_mea;
		
		/* --- Hatta number and enhancement factor --- */
        Ha = pow((k_r * d_co2), 0.5) / kL;              /* Hatta number          [-] */
        E  = pow((1.0 + pow(Ha, 2)), 0.5);               /* enhancement factor (fast reaction approximation) */
		
		/* --- Henry's constant --- */
        Henry = exp(896.5 / T_l - 10.5913);              /* [mol/(m³·kPa)] */
		
		/* --- driving force yL = equilibrium CO2 conc. - actual CO2 conc.  [mol/m³] --- */
        yL = exp(35.81 * y_meacoo - 3.09) * Henry
             - (C_P(cell, thread_gas) / 1000.0 + 101.325) * Henry
             * (C_YI(cell, thread_gas, 0)
             / (C_YI(cell, thread_gas, 0) + 2.44 * (1.0 - C_YI(cell, thread_gas, 0))));
			 
		 /* --- MEAH liq-phase mass source term --- */
        /* -1 * E * kL [m/s] * ae [m²/m³] * M_MEAH [kg/mol] * yL [mol/m³] = [kg/(m³·s)] */
        m_dot_meah = - E * kL * ae * 62.0 * yL;
		
		   dS[eqn] = 0.0;                
    }
    else
    {
        /* source set to zero when VOF is outside the valid interface range */
        m_dot_meah = 0.0;
        dS[eqn]   = 0.0;
    }

    return m_dot_meah;
}	


/*------------------------------------------------------------
 * DEFINE_SOURCE : meacoo_src
 * Returns : meacoo mass source term in the Liq phase  [kg/(m³·s)]
 *           Positive -> liquid generates MEACOO (absorption)
 *           Negative -> liquid consumes MEACOO (desorption)
 *------------------------------------------------------------*/
 DEFINE_SOURCE(liqmeacoo_src, cell, thread_liq, dS, eqn)
{
    Thread *thread_mix, *thread_gas;
    thread_mix = THREAD_SUPER_THREAD(thread_liq);        /* mixture thread              */
    thread_gas = THREAD_SUB_THREAD(thread_mix, 1);       /* gas sub-thread (index 1) */
 
    real ae, c_mea, d_n2o, d_co2, c_d, kL, k_r, Ha, E;
    real yL, m_dot_meacoo, Henry, loading, v, c_meacoo;
 
    real T_l = C_T(cell, thread_liq);                    /* liquid temperature  [K] */
	if (C_VOF(cell, thread_gas) < 0.99 && C_VOF(cell, thread_gas) > 0.01)
    {
        /* --- interfacial area density (spherical bubble assumption) --- */
        /* ae = 6 * alpha_g / d_b  [m²/m³] */
        ae = 6.0 * C_VOF(cell, thread_gas) / C_PHASE_DIAMETER(cell, thread_gas);
		
		 /* --- species and CO2 loading --- */
        real y_meacoo    = C_YI(cell, thread_liq, 2);          /* MEACOO- mass fraction  [-]     */
        real denom = 104.08 - 44.01 * y_meacoo;
        loading    = y_meacoo * 203.6 / denom;                 /* CO2 loading  [mol/mol]         */
		
		/* --- liquid dynamic viscosity --- */
        v = exp(7.08 - 0.021 * T_l + 0.97 * loading) / 1000.0;     /* [Pa·s] */
		
		/* --- molar concentrations --- */
        c_mea    = C_YI(cell, thread_liq, 0) * density / 61.0;      /* MEA     [mol/L] */
        c_meacoo = C_YI(cell, thread_liq, 2) * density / 104.0;     /* MEACOO- [mol/L] */
		
		/* --- diffusivity via N2O analogy --- */
        d_n2o = (5.07 + 0.865 * c_mea + 0.278 * c_mea * c_mea)
                * exp((-2371.0 - 93.4 * c_mea) / T_l) * 1e-6;      /* D_N2O  [m²/s] */
        d_co2 = d_n2o * 0.4635 * exp(252.0 / T_l);                  /* D_CO2  [m²/s] */
		
		/* --- turbulent energy dissipation rate --- */
        /* epsilon = C_mu * k * omega  [m²/s³] */
        c_d = 0.09 * C_K(cell, thread_liq) * C_O(cell, thread_liq);
		
		/* --- liquid-side mass transfer coefficient (Calderbank-Moo-Young) --- */
        kL = K * pow(d_co2, 0.5) * pow((c_d / v * density), 0.25); /* [m/s] */
		
		 /* --- pseudo-first-order reaction rate constant --- */
        /* k_r = k2 * c_MEA  [1/s] */
        k_r = pow(10.0, (10.99 - 2152.0 / T_l)) * c_mea;
		
		/* --- Hatta number and enhancement factor --- */
        Ha = pow((k_r * d_co2), 0.5) / kL;              /* Hatta number          [-] */
        E  = pow((1.0 + pow(Ha, 2)), 0.5);               /* enhancement factor (fast reaction approximation) */
		
		/* --- Henry's constant --- */
        Henry = exp(896.5 / T_l - 10.5913);              /* [mol/(m³·kPa)] */
		
		/* --- driving force yL = equilibrium CO2 conc. - actual CO2 conc.  [mol/m³] --- */
        yL = exp(35.81 * y_meacoo - 3.09) * Henry
             - (C_P(cell, thread_gas) / 1000.0 + 101.325) * Henry
             * (C_YI(cell, thread_gas, 0)
             / (C_YI(cell, thread_gas, 0) + 2.44 * (1.0 - C_YI(cell, thread_gas, 0))));
			 
		 /* --- MEACOO liq-phase mass source term --- */
        /* -1 * E * kL [m/s] * ae [m²/m³] * M_MEACOO [kg/mol] * yL [mol/m³] = [kg/(m³·s)] */
        m_dot_meacoo = - E * kL * ae * 104.0 * yL;
		
		   dS[eqn] = 0.0;                
    }
    else
    {
        /* source set to zero when VOF is outside the valid interface range */
        m_dot_meacoo = 0.0;
        dS[eqn]   = 0.0;
    }

    return m_dot_meacoo;
}	


/*------------------------------------------------------------
 * DEFINE_SOURCE : liqmea_heat_src
 * Returns : reaction heat source in the LIQUID phase  [W/m³]
 * Note    : Reads the MEA mass source from UDMI[2] to back-calculate
 *           the CO2 molar flux, then multiplies by the reaction enthalpy.
 *           Negative sign: desorption is endothermic (absorption is exothermic).
 *------------------------------------------------------------*/
DEFINE_SOURCE(liqmea_heat_src, cell, thread_liq, dS, eqn)
{
    real source = 0.0;
    real m_dot_co2 = C_UDMI(cell, thread_liq, 2);       /* CO2 mass source  [kg/(m³·s)] from UDMI[2] */
	
	real M_CO2 = 0.044;                                  /* molar mass used in this term [kg/mol]*/
	
	real reaction_rate_co2_mol = m_dot_co2 / M_CO2;          /* [mol/(m³·s)] */
	
	real heat_of_reaction_j_per_mol = 75000.0;           /* reaction enthalpy ~75 kJ/mol for MEA-CO2 */
	
	    /* negative: heat is consumed during desorption */
    source = -reaction_rate_co2_mol * heat_of_reaction_j_per_mol;  /* [W/m³] */
 
    dS[eqn] = 0.0;
    return source;
}


/*------------------------------------------------------------
 * DEFINE_SOURCE : liq_src
 * Returns : total liquid-phase mass source term  [kg/(m³·s)]
 * Note    : Accounts for liquid mass lost due to CO2 leaving the liquid.
 *           M_CO2 = 44 g/mol; negative sign = liquid phase loses mass.
 *------------------------------------------------------------*/
 
 DEFINE_SOURCE(liq_src, cell, thread_liq, dS, eqn)
{
    Thread *thread_mix, *thread_gas;
    thread_mix = THREAD_SUPER_THREAD(thread_liq);        /* mixture thread              */
    thread_gas = THREAD_SUB_THREAD(thread_mix, 1);       /* gas sub-thread (index 1) */
 
    real ae, c_mea, d_n2o, d_co2, c_d, kL, k_r, Ha, E;
    real yL, m_dot_liq, Henry, loading, v, c_meacoo;
 
    real T_l = C_T(cell, thread_liq);                    /* liquid temperature  [K] */
	if (C_VOF(cell, thread_gas) < 0.99 && C_VOF(cell, thread_gas) > 0.01)
    {
        /* --- interfacial area density (spherical bubble assumption) --- */
        /* ae = 6 * alpha_g / d_b  [m²/m³] */
        ae = 6.0 * C_VOF(cell, thread_gas) / C_PHASE_DIAMETER(cell, thread_gas);
		
		 /* --- species and CO2 loading --- */
        real y_meacoo    = C_YI(cell, thread_liq, 2);          /* MEACOO- mass fraction  [-]     */
        real denom = 104.08 - 44.01 * y_meacoo;
        loading    = y_meacoo * 203.6 / denom;                 /* CO2 loading  [mol/mol]         */
		
		/* --- liquid dynamic viscosity --- */
        v = exp(7.08 - 0.021 * T_l + 0.97 * loading) / 1000.0;     /* [Pa·s] */
		
		/* --- molar concentrations --- */
        c_mea    = C_YI(cell, thread_liq, 0) * density / 61.0;      /* MEA     [mol/L] */
        c_meacoo = C_YI(cell, thread_liq, 2) * density / 104.0;     /* MEACOO- [mol/L] */
		
		/* --- diffusivity via N2O analogy --- */
        d_n2o = (5.07 + 0.865 * c_mea + 0.278 * c_mea * c_mea)
                * exp((-2371.0 - 93.4 * c_mea) / T_l) * 1e-6;      /* D_N2O  [m²/s] */
        d_co2 = d_n2o * 0.4635 * exp(252.0 / T_l);                  /* D_CO2  [m²/s] */
		
		/* --- turbulent energy dissipation rate --- */
        /* epsilon = C_mu * k * omega  [m²/s³] */
        c_d = 0.09 * C_K(cell, thread_liq) * C_O(cell, thread_liq);
		
		/* --- liquid-side mass transfer coefficient (Calderbank-Moo-Young) --- */
        kL = K * pow(d_co2, 0.5) * pow((c_d / v * density), 0.25); /* [m/s] */
		
		 /* --- pseudo-first-order reaction rate constant --- */
        /* k_r = k2 * c_MEA  [1/s] */
        k_r = pow(10.0, (10.99 - 2152.0 / T_l)) * c_mea;
		
		/* --- Hatta number and enhancement factor --- */
        Ha = pow((k_r * d_co2), 0.5) / kL;              /* Hatta number          [-] */
        E  = pow((1.0 + pow(Ha, 2)), 0.5);               /* enhancement factor (fast reaction approximation) */
		
		/* --- Henry's constant --- */
        Henry = exp(896.5 / T_l - 10.5913);              /* [mol/(m³·kPa)] */
		
		/* --- driving force yL = equilibrium CO2 conc. - actual CO2 conc.  [mol/m³] --- */
        yL = exp(35.81 * y_meacoo - 3.09) * Henry
             - (C_P(cell, thread_gas) / 1000.0 + 101.325) * Henry
             * (C_YI(cell, thread_gas, 0)
             / (C_YI(cell, thread_gas, 0) + 2.44 * (1.0 - C_YI(cell, thread_gas, 0))));
			 
		/* only CO2 leaves the liquid phase; M_CO2 = 44 g/mol */
        m_dot_liq = - E * kL * ae * 44 * yL;
		
		   dS[eqn] = 0.0;                
    }
    else
    {
        /* source set to zero when VOF is outside the valid interface range */
        m_dot_liq = 0.0;
        dS[eqn]   = 0.0;
    }

    return m_dot_liq;
}	

