#parameters
s2o.par <- list(
  am = list(ic = c('OO'=2.1, 'PP'=1.9),
            pp = c("p0"=1.0, "d0"=1.0, "p3"=1.0, "d3"=1.0, "bp0"=0.05, "bd0"=0.05, "bp3"=0.05, "bd3"=0.05, "nt"=1, "nd"=2, "tot"=4)),
  
  ti = list(ic = c('OO'=2.1, 'PO'=0, 'PP'=1.9),
            pp = c("p0"=1.0, "p1"=1.0, "p2"=1.0, "p3"=1.0, 
                   "d0"=1.0, "d1"=1.0, "d2"=1.0, "d3"=1.0,
                   "bp0"=0.05, "bp1"=0.05, "bp2"=0.05, "bp3"=0.05,
                   "bd0"=0.05, "bd1"=0.05, "bd2"=0.05, "bd3"=0.05,
                   "nt"=1, "nd"=2, "tot"=4)),
  
  ai = list(ic = c('OO'=2.1, 'PO'=0, 'PP'=1.9),
            pp = c("p0"=1.0, "p1"=1.0, "p2"=0.0, "p3"=1.0, 
                   "d0"=1.0, "d1"=0.0, "d2"=1.0, "d3"=1.0, 
                   "bp0"=0.05, "bp1"=0.05, "bp2"=0.0, "bp3"=0.05,
                   "bd0"=0.05, "bd1"=0.0,"bd2"=0.05, "bd3"=0.05, 
                   "nt"=1, "nd"=2, "tot"=4)),
  
  bp = list(ic = c('OO'=2.1, 'PO'=0, 'PP'=1.9),
            pp = c("p0"=1.0, "p1"=1.0, "p2"=0.0, "p3"=1.0, 
                   "d0"=0.0, "d1"=1.0, "d2"=1.0, "d3"=1.0,
                   "bp0"=0.05, "bp1"=0.05, "bp2"=0.0, "bp3"=0.05,
                   "bd0"=0.0, "bd1"=0.05, "bd2"=0.05, "bd3"=0.5,
                   "nt"=1, "nd"=2, "tot"=4)),
  
  bd1 = list(ic = c('OO'=2.1, 'PO'=0, 'PP'=1.9),
             pp = c("p0"=1.0, "p1"=1.0, "p2"=0.0, "p3"=1.0, 
                    "d0"=0.0, "d1"=1.0, "d2"=1.0, "d3"=1.0,
                    "bp0"=0.05, "bp1"=0.05, "bp2"=0.0, "bp3"=0.05,
                    "bd0"=0.0, "bd1"=0.0, "bd2"=0.05, "bd3"=0.5,
                    "nt"=1, "nd"=2, "tot"=4)),
  
  so = list(ic = c('OO'=2.1, 'PO'=0, 'PP'=1.9),
            pp = c("p0"=1.0, "p1"=0.0, "p2"=0.0, "p3"=0.0, 
                   "d0"=0.0, "d1"=1.0, "d2"=0.0, "d3"=0.0,
                   "bp0"=0.0, "bp1"=0.0, "bp2"=0.0, "bp3"=0.05,
                   "bd0"=0.0, "bd1"=0.0, "bd2"=0.05, "bd3"=0.0,
                   "nt"=1, "nd"=2, "tot"=4)),
  
  co = list(ic = c('OO'=2.1, 'PO'=0, 'PP'=1.9),
            pp = c("p0"=1.0, "p1"=0.0, "p2"=0.0, "p3"=1.0, 
                   "d0"=0.0, "d1"=1.0, "d2"=0.0, "d3"=0.0,
                   "bp0"=0.00, "bp1"=0.0, "bp2"=0.0, "bp3"=0.0,
                   "bd0"=0.0, "bd1"=0.0, "bd2"=0.05, "bd3"=0.0,
                   "nt"=1, "nd"=2, "tot"=4))
  
)#end list




#models
s2o.mod <- list(
  
  f1 = function(tt,x,k){
    with(as.list(c(x,k)),{
      dOO = -p0*nt*PP*OO + d0*nd*OO*(tot-PP-OO) - bp0*nt*OO + bd0*nd*(tot-PP-OO)
      dPP = p3*nt*PP*(tot-PP-OO) - d3*nd*OO*PP + bp3*nt*(tot-PP-OO) - bd3*nd*PP
      return(list(c(dOO,dPP)))
    })
  },
  
  f5 = function(tt, x, k){
    with(as.list(c(x, k)),{
      dOO = -p0*nt*PP*OO - p1*nt*PP*OO + d0*nd*OO*(tot-PP-OO-PO) + d1*nd*OO*PO - bp0*nt*OO - bp1*nt*OO + bd0*nd*(tot-PP-OO-PO) + bd1*nd*PO
      dPP = p3*nt*PP*(tot-PP-OO-PO) + p2*nt*PP*PO - d2*nd*OO*PP - d3*nd*OO*PP + bp3*nt*(tot-PP-OO-PO) + bp2*nt*PO - bd2*nd*PP - bd3*nd*PP
      dPO = p1*nt*PP*OO - p2*nt*PP*PO + d2*nd*OO*PP - d1*nd*OO*PO + bp1*nt*OO - bp2*nt*PO + bd2*nd*PP - bd1*nd*PO
      return (list (c(dOO, dPO, dPP)))
    })}
  
)#end list
