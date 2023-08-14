# Import Packages and Libraries
import numpy as np
import math
import pandas as pd

# Define the LewisFun function for tray calculations
def LewisFun(F, D, r, xf, xd, xb, xq, Tbe, Tbw, Fm, Dm, Bm, Tbmixf):  
    # Parameters (constants)
    Tray_eff = 0.7  # efficiency of trays
    TS = 0.6  # distance between trays [m]
    rhoe = 789  # ethanol density [kg/m3]
    rhow = 1000  # water density [kg/m3]
    Cpe = 2.46  # heat capacity of ethanol [kJ/kg oC]
    Cpw = 4.186  # heat capacity of water [kJ/kg oC]
    He = 855  # latent heat of ethanol [kJ/kg]
    Hw = 2264.705  # latent heat of water [kJ/kg]
    Hs = 2080  # latent heat of steam [kJ/kg]
    Twout = 30  # cooling water outlet temperature [oC]
    Twin = 20  # cooling water inlet temperature [oC]
    Scost = 0.03  # cost per kg steam [$]
    Wcost = 1.48*10**(-5)  # cost per kg cooling water [$]
    CEPCI68 = 114  # Chemical Engineering Plant Cost Index at 1968
    CEPCI24 = 715  # Chemical Engineering Plant Cost Index at 2024
    IR = 1.0328  # inflation rate average in US
    Time = 25  # life time of process unit, [yr]
    
    # Lewis Structure
    x = np.array([])  # 1D array for the liquid ethanol molar fraction
    y = np.array([])  # 1D array for the vapor ethanol molar fraction

    # Rectifying section calculations
    n = 0  # indexing
    y = np.append(y, xd)  # initial condition
    Tbmix = xd*Tbe + (1-xd)*Tbw  # boiling point of distillate [oC]
    a = 0.6078*Tbmix-46.373  # relative volatility derived via linear regression of experimental data
    x_element = y[n]/(a-(a-1)*y[n])  # vapour liquid equilibrium
    x = np.append(x, x_element)
    while (x[n]>xq):
        y_element = r/(r+1)*x[n] + 1/(r+1)*xd
        y = np.append(y, y_element)
        Tbmix = x[n]*Tbe + (1-x[n])*Tbw
        a = 0.6078*Tbmix-46.373
        x_element = y[n+1]/(a-(a-1)*y[n+1])
        x = np.append(x, x_element)
        n = n+1    
    Nr = y.size-1  # number of theoretical trays in rectifying section

    # Repetitive Structure for calculations in stripping section
    s = (xb-xf)/(xf-xd)*(r+1)  # Boil-up ratio
    while (x[n]>xb):
        y_element = (s+1)/s*x[n] - 1/s*xb
        y = np.append(y, y_element)
        Tbmix = x[n]*Tbe + (1-x[n])*Tbw
        a = 0.6078*Tbmix-46.373
        x_element = y[n+1]/(a-(a-1)*y[n+1])
        x = np.append(x, x_element)
        n = n+1
    Ns = y.size-2 + (x[n-1]-xb)/(x[n-1]-x[n]) - Nr  # number of theoretical trays in stripping section

    # Mixture Density in Bottoms
    rhomix = 1/((xb/rhoe)+((1-xb)/rhow))  # density of mixture in bottoms

    # Calculation of Column's Characteristics
    N = y.size-1 + (x[n-1]-xb)/(x[n-1]-x[n])  # number of theoretical trays
    N_real = np.ceil(N/Tray_eff)  # number of trays
    Nr_real = np.ceil(Nr/Tray_eff)  # number of trays in rectifying section
    Ns_real = np.ceil(Ns/Tray_eff)  # number of trays in stripping section
    Fp = Nr_real + 0.5  # feed position
    Hc = 1.2*TS*(N_real-1)  # height of column [m]
    V = D+r*D  # vapor molar flow rate in column [mol/h]
    Dc = (4*V/(0.6*math.pi*rhomix))**(1/2)  # diameter of column [m]

    # Energy Balances Calculations
    Hd = He*xd + Hw*(1-xd)  # latent heat of distillate [kJ/kg]
    Hb = He*xb + Hw*(1-xb)  # latent heat of bottoms [kJ/kg]
    Cpd = Cpe*xd + Cpw*(1-xd)  # heat capacity of distillate [kJ/kg oC]
    Cpf = Cpe*xf + Cpw*(1-xf)  # heat capacity of distillate [kJ/kg oC]
    Tbmixb = xb*Tbe + (1-xb)*Tbw
    Qr = Fm*Cpf*(Tbmixb-Tbmixf) + V*Hd  # reboiler's duty [kJ/h]
    Tbmixd = xd*Tbe + (1-xd)*Tbw
    Qc = V*Hd  # condenser's duty [kJ/h]
    ms = Qr/Hs  # mass of steam needed [kg/h]
    mw = Qc/(Cpw*(Twout-Twin))  # mass of cooling water needed [kg/h]

    # Calculation of Cost
    CBMC = (CEPCI24/CEPCI68)*5.5*935.6*Hc**0.81*Dc**1.05  # cost of column's shell [$]
    CBMT = (CEPCI24/CEPCI68)*1.7*125.2*(Hc/1.2)**0.97*Dc**1.45  # cost of column's trays [$]
    UCS = Scost*ms*365*24*0.95  # annual cost of steam [$/yr]
    UCW = Wcost*mw*365*24*0.95  # annual cost of cooling water [$/yr]
    CEB = 1.28*1000/rhoe*Bm*xb*7920  # annual revenue loss of ethanol in bottoms [$/yr]
    TOP = UCS+UCW+CEB  # total operating cost [$/plant lifetime]
    t = 1
    while t<Time:
        TOP = TOP + (UCS+UCW+CEB)*IR
        t = t+1
    TC = CBMC+CBMT+TOP  # total cost of distillation column [$/plant lifetime]
   
    return s, N_real, Nr_real, Ns_real, Fp, Dc, Hc, Tbmixd, Tbmixb, Qc, Qr, ms, mw, CBMC, CBMT, UCS, UCW, CEB, TC

# Define the ColumnDesignFun function for column design calculations
def ColumnDesignFun(Dme, xfm, xdm, xbm):
    # Parameters (constants)
    q = 1  # parameter equals 1 if feed is saturated liquid
    MWe = 18.015  # molecular mass of ethanol [kg/kmol]
    MWw = 46.07  # molecular mass of water [kg/kmol]
    Tbe = 78.37  # boiling point of ethanol [oC]
    Tbw = 100  # boiling point of water [oC]
    
    Dm = Dme/xdm  # convert produced ethanol mass flow rate into total distillate mass flow rate [kg/h]
    
    # Mass to mole conversion
    Med = Dm*xdm
    Ned = Med/MWe
    Mwd = Dm*(1-xdm)
    Nwd = Mwd/MWw
    D = Ned+Nwd
    xf = (xfm/MWe)/((xfm/MWe)+((1-xfm)/MWw))
    xq = xf*q
    xb = (xbm/MWe)/((xbm/MWe)+((1-xbm)/MWw))
    xd = (xdm/MWe)/((xdm/MWe)+((1-xdm)/MWw))
    
    # Mass Balances Calculations
    F = D*(xd-xb)/(xf-xb)  # feed molar flow rate [kmol/h]
    B = F - D  # bottoms molar flow rate [kmol/h]
    De = D*xd  # ethanol in distillate [mol]
    Be = B*xb  # ethanol in bottoms [mol]
    Dw = D*(1-xd)  # ethanol in distillate [mol]
    Bw = B*(1-xb)  # water in bottoms [mol]
    DRe = De/Be  # ethanol distribution ration
    DRw = Dw/Bw  # water distribution ration
    
    # Mass to mole conversion
    Neb = B*xb
    Meb = Neb*MWe
    Nwb = B*(1-xb)
    Mwb = Nwb*MWw
    Bm = Meb+Mwb
    Nef = B*xf
    Mef = Nef*MWe
    Nwf = B*(1-xf)
    Mwf = Nwf*MWw
    Fm = Mef+Mwf

    # Average Relative Volatility Calculations
    Tbmixf = xf*Tbe + (1-xf)*Tbw
    a0 = 0.6078*Tbmixf-46.373

    # Minimum Reflux & Number of Trays Calculations
    Nmin = np.log10(xd/(1-xd)*(1-xb)/xb)/np.log10(a0)  # Fenske minimum number of trays equation
    if q == 1:
        theta = a0/(1-xf+a0*xf)  # Underwood theta equation
    rmin = (a0*xd)/(a0-theta)+(1-xd)/(1-theta)-1  # Underwood minimum reflux rate equation
    
    # Find the Optimal Reflux Ratio, Number of Trays and Feed Position for known Mass Fractions and Demand
    columns = ['Feed mass flow (kg/h)', 'Distillate mass flow (kg/h)', 'Bottoms mass flow (kg/h)',
               'Feed mass fraction', 'Distillate mass fraction', 'Bottoms mass fraction',  
               'R/Rmin', 'Reflux Ratio', 'Boil-up Ratio', 'Number of Trays', 'Number of Trays in Rectifying Section', 
               'Number of Trays in Stripping Section', 'Feed location', 'Diameter (m)', 'Height (m)', 
               'Condenser Temperature (oC)', 'Reboiler Temperature (oC)', 'Condenser Duty (kJ/h)', 
               'Reboiler duty (kJ/h)', 'Steam mass (kg/h)', 'Cooling Water mass (kg/h)', 'CBMC ($)', 'CBMT ($)', 
               'UCS ($)', 'UCW ($)', 'CEB ($)', 'TC ($)']
    df_r = pd.DataFrame(columns=columns)  # create an empty dataframe with specified columns
    
    if rmin <= 0:  # rmin positive check
        print("      Error message: Infeasible Column")
        new_row = {'Feed mass flow (kg/h)': math.nan, 'Distillate mass flow (kg/h)': math.nan, 'Bottoms mass flow (kg/h)': math.nan,
                   'Feed mass fraction': math.nan, 'Distillate mass fraction': math.nan, 'Bottoms mass fraction': math.nan,  
                   'R/Rmin': math.nan, 'Reflux Ratio': math.nan, 'Boil-up Ratio': math.nan, 'Number of Trays': math.nan, 
                   'Number of Trays in Rectifying Section': math.nan, 'Number of Trays in Stripping Section': math.nan, 
                   'Feed location': math.nan, 'Diameter (m)': math.nan, 'Height (m)': math.nan, 'Condenser Temperature (oC)': math.nan, 
                   'Reboiler Temperature (oC)': math.nan, 'Condenser Duty (kJ/h)': math.nan, 'Reboiler duty (kJ/h)': math.nan, 
                   'Steam mass (kg/h)': math.nan, 'Cooling Water mass (kg/h)': math.nan, 'CBMC ($)': math.nan, 'CBMT ($)': math.nan, 
                   'UCS ($)': math.nan, 'UCW ($)': math.nan, 'CEB ($)': math.nan, 'TC ($)': math.nan}
        min_row = pd.DataFrame([new_row], index=[0])
    else:
        r = 1.01*rmin  # initial r value
        while (r<=3*rmin):  # range of r values
            rr = r/rmin

            s, N_real, Nr_real, Ns_real, Fp, Dc, Hc, Tbmixd, Tbmixb, Qc, Qr, ms, mw, CBMC, CBMT, UCS, UCW, CEB, TC = LewisFun(F, D, r, xf, xd, xb, xq, Tbe, Tbw, Fm, Dm, Bm, Tbmixf)

            # Add Results to Dataframe df_r
            new_row = {'Feed mass flow (kg/h)': Fm, 'Distillate mass flow (kg/h)': Dm, 'Bottoms mass flow (kg/h)': Bm,
                   'Feed mass fraction': xfm, 'Distillate mass fraction': xdm, 'Bottoms mass fraction': xbm,  
                   'R/Rmin': rr, 'Reflux Ratio': r, 'Boil-up Ratio': s, 'Number of Trays': N_real, 
                   'Number of Trays in Rectifying Section': Nr_real, 'Number of Trays in Stripping Section': Ns_real, 
                   'Feed location': Fp, 'Diameter (m)': Dc, 'Height (m)': Hc, 'Condenser Temperature (oC)': Tbmixd, 
                   'Reboiler Temperature (oC)': Tbmixb, 'Condenser Duty (kJ/h)': Qc, 'Reboiler duty (kJ/h)': Qr, 
                   'Steam mass (kg/h)': ms, 'Cooling Water mass (kg/h)': mw, 'CBMC ($)': CBMC, 'CBMT ($)': CBMT, 
                   'UCS ($)': UCS, 'UCW ($)': UCW, 'CEB ($)': CEB, 'TC ($)': TC}
            df_r = pd.concat([pd.DataFrame(new_row, index=[0]), df_r], ignore_index=True)  # add a new row with the specified values
            r = r+0.01*rmin
        
        # Find the Column with the Minimum Total Annual Cost
        min_value = df_r['TC ($)'].min()
        min_row = df_r.loc[df_r['TC ($)'] == min_value]  # find the row that has the minimum value

    return min_row


# Parameter
Dme = 33871  # ethanol produced mass flow [kg/h]
xfm = 0.01  # feed volatile compound mass fraction, initialization
xdm = 0.4  # distillate volatile compound mass fraction, initialization
xbm = 0.001  # bottoms volatile compound mass fraction, initialization

# Find the Optimal Reflux Ratio, Number of Trays and Feed Position for known Mass Fractions and Demand
columns = ['Feed mass flow (kg/h)', 'Distillate mass flow (kg/h)', 'Bottoms mass flow (kg/h)',
               'Feed mass fraction', 'Distillate mass fraction', 'Bottoms mass fraction',  
               'R/Rmin', 'Reflux Ratio', 'Boil-up Ratio', 'Number of Trays', 'Number of Trays in Rectifying Section', 
               'Number of Trays in Stripping Section', 'Feed location', 'Diameter (m)', 'Height (m)', 
               'Condenser Temperature (oC)', 'Reboiler Temperature (oC)', 'Condenser Duty (kJ/h)', 
               'Reboiler duty (kJ/h)', 'Steam mass (kg/h)', 'Cooling Water mass (kg/h)', 'CBMC ($)', 'CBMT ($)', 
               'UCS ($)', 'UCW ($)', 'CEB ($)', 'TC ($)']
df_c = pd.DataFrame(columns=columns)  # create an empty dataframe with specified columns
    
while (xfm<=0.2):  # range of xfm values
    print(f'Execute outer loop, xfm = {xfm}')
    while (xdm<=0.9):  # range of xdm values
        print(f'    Execute middle loop, xdm = {xdm}')
        while (xbm<=0.01):  # range of xbm values
            print(f'        Execute inner loop, xbm = {xbm}')
            min_row = ColumnDesignFun(Dme, xfm, xdm, xbm)
            min_row.reset_index(drop=True, inplace=True)  # reset index to match df_c
            df_c = pd.concat([df_c, min_row], ignore_index=True)
            xbm = xbm + 0.001
        xdm = xdm + 0.01
        xbm = 0.001  # bottoms mass fraction initialization
        print('------------------------------------------------------------------')
    xfm = xfm + 0.01
    xdm = 0.4  # distillate mass fraction initialization 
    print('------------------------------------------------------------------')

df_c = df_c.dropna()
df_c.to_excel('distcol.xlsx', index=False)
print("Distillation Column Converged Normally --> all data in distcol.xlsx file")
