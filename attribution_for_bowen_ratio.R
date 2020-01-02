rm(list = ls())
# loading required libraries 
library(Deriv)

# Please refer a paper 000 for the details 

## Delta Bowen ratio
# partial derivative for Bowen Ratio
bowen <- function(Ta,Press,SWin,LWin,ra,rs,emis,qa,albedo,G){
  rho_air = Press/(r_d*Ta) # density of air
  Ta_c = Ta-273.15
  lambda_o = 1/(4*emis*sb*Ta^3)
  Rn_sat = SWin*(1-albedo)+emis*LWin-emis*sb*Ta^4
  qa_sat = (0.622/Press)*611*exp((17.27*Ta_c)/(Ta_c+237.3))
  diff_qasat_Ta=(0.622/Press)*exp((17.27*Ta_c)/(Ta_c+237.3))*(2508.3/(Ta_c+237.3)^2)*1000
  Ts = ((Rn_sat-G)-rho_air*lv/(ra+rs)*(qa_sat-qa))/((1/lambda_o)+rho_air*cp/ra+rho_air*lv/(ra+rs)*diff_qasat_Ta)+Ta
  
  H = rho_air*cp/ra*(Ts-Ta)
  LE = rho_air*lv/(ra+rs)*(qa_sat+diff_qasat_Ta*(Ts-Ta)-qa)
  H/LE
}
diff_Ta     <- Deriv(bowen,'Ta')
diff_Press  <- Deriv(bowen,'Press')
diff_SWin   <- Deriv(bowen,'SWin')
diff_LWin   <- Deriv(bowen,'LWin')
diff_ra     <- Deriv(bowen,'ra')
diff_rs     <- Deriv(bowen,'rs')
diff_qa     <- Deriv(bowen,'qa')
diff_albedo <- Deriv(bowen,'albedo')
diff_G      <- Deriv(bowen,'G')

# calculations
dBn_dswd  = diff_SWin(Ta,Press,SWin,LWin,ra,rs,emis,qa,albedo,G)
dBn_dalpha = diff_albedo(Ta,Press,SWin,LWin,ra,rs,emis,qa,albedo,G)
dBn_drld = diff_LWin(Ta,Press,SWin,LWin,ra,rs,emis,qa,albedo,G)
dBn_dra = diff_ra(Ta,Press,SWin,LWin,ra,rs,emis,qa,albedo,G)
dBn_drs = diff_rs(Ta,Press,SWin,LWin,ra,rs,emis,qa,albedo,G)
dBn_dqa = diff_qa(Ta,Press,SWin,LWin,ra,rs,emis,qa,albedo,G)
dBn_dTa = diff_Ta(Ta,Press,SWin,LWin,ra,rs,emis,qa,albedo,G)
dBn_dGrnd = diff_G(Ta,Press,SWin,LWin,ra,rs,emis,qa,albedo,G)


