library(EpiModel)
# program for default model, no change after finalized

# customized SEIR with quarantine and imported cases;

SEIRQ <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # aggregate parameters over young and old
    num = ynum + snum
    q.num = qy.num + qs.num
    s.num = sy.num + ss.num
    e.num = ey.num + es.num
    #total reported prev cases;
    i.num = iyd.num + iyq.num + isd.num + isq.num
    id.num = iyd.num + isd.num
    iu.num = iyu.num + isu.num
    iq.num = iyq.num + isq.num
    h.num = hy.num + hs.num
    d.num = dy.num + ds.num
    r.num =ry.num + rs.num
    
    # premature stop importing cases at 20 days
    if(t >timpt) fy2 = 0
    
    dNY <- (fy1 + fy2)*ynum  
    dNS <- (fs1 + fs2)*snum 
    dQY <- qy*(i.num/s.num)*sy.num - (p + gy + uy)*qy.num
    dQS <- qs*(i.num/s.num)*ss.num - (p + gs + us)*qs.num
    dSY <- p*qy.num + fy1*ynum - ((by1*iyd.num + by2*iyu.num)*ky + (bs1*isd.num + bs2*isu.num)*kys)*sy.num/ynum - (qy*(i.num/s.num) + uy)*sy.num
    dSS <- p*qs.num + fs1*snum - ((bs1*isd.num + bs2*isu.num)*ks + (by1*iyd.num + by2*iyu.num)*kys)*ss.num/snum - (qs*(i.num/s.num) + us)*ss.num
    dEY <- fy2*ynum + ((by1*iyd.num + by2*iyu.num)*ky + (bs1*isd.num + bs2*isu.num)*kys)*sy.num/ynum - (vy1 + vy2 + uy)*ey.num
    dES <- fs2*snum + ((bs1*isd.num + bs2*isu.num)*ks + (by1*iyd.num + by2*iyu.num)*kys)*ss.num/snum - (vs1 + vs2 + us)*es.num
    dIYQ <- gy*qy.num - (hy1 +ry1)*iyq.num
    dISQ <- gs*qs.num - (hs1 +rs1)*isq.num
    dIYD <- vy1*ey.num  - (hy2 + ry2)*iyd.num
    dISD <- vs1*es.num  - (hs2 + rs2)*isd.num
    dIYU <- vy2*ey.num  - (ry3 + uy)*iyu.num
    dISU <- vs2*es.num  - (rs3 + us)*isu.num
    dHY <- hy1*iyq.num + hy2*iyd.num -(ry4 + wy)*hy.num
    dHS <- hs1*isq.num + hs2*isd.num -(rs4 + ws)*hs.num
    dRY <- ry1*iyq.num + ry2*iyd.num + ry3*iyu.num + ry4*hy.num - uy*ry.num
    dRS <- rs1*isq.num + rs2*isd.num + rs3*isu.num + rs4*hs.num - us*rs.num
    dDY <- wy*hy.num
    dDS <- ws*hs.num
    
    
    # Compartments and flows are part of the derivative vector, must be in the same order as para and control
    # Other calculations to be output are outside the vector, but within the containing list
    list(c(dNY,dNS,dQY,dQS,dSY,dSS,dEY,dES,dIYQ,dISQ,dIYD,dISD,dIYU,dISU,dHY,dHS, dRY,dRS,dDY,dDS,
           ns.flow =  fy1*ynum + fs1*snum,
           sqy.flow = qy*(i.num/s.num)*sy.num,
           sqs.flow = qs*(i.num/s.num)*ss.num,
           qsy.flow = p*qy.num,
           qss.flow = p*qs.num,
           sey.flow = fy2*ynum + ((by1*iyd.num + by2*iyu.num)*ky + (bs1*isd.num + bs2*isu.num)*kys)*sy.num/ynum,
           ses.flow = fs2*snum + ((bs1*isd.num + bs2*isu.num)*ks + (by1*iyd.num + by2*iyu.num)*kys)*ss.num/snum,
           ney.flow = fy2*ynum,
           nes.flow = fs2*snum,
           eiyd.flow = vy1*ey.num,
           eisd.flow = vs1*es.num,
           eiyu.flow = vy2*ey.num,
           eisu.flow = vs2*es.num,
           qiyq.flow = gy*qy.num,
           qisq.flow = gs*qs.num,
           iyqh.flow = hy1*iyq.num,
           isqh.flow = hs1*isq.num,
           iyqr.flow = ry1*iyq.num,
           isqr.flow = rs1*isq.num,
           iydh.flow = hy2*iyd.num,
           isdh.flow = hs2*isd.num,
           iydr.flow = ry2*iyd.num,
           isdr.flow = rs2*isd.num,
           iyur.flow = ry3*iyu.num,
           isur.flow = rs3*isu.num,
           hyr.flow = ry4*hy.num,
           hsr.flow = rs4*hs.num,
           hyd.flow = wy*hy.num,
           hsd.flow = ws*hs.num,
           hy.flow = hy1*iyq.num +hy2*iyd.num,
           hs.flow = hs1*isq.num +hs2*isd.num,
           iy.flow = vy1*ey.num + gy*qy.num,
           is.flow = vs1*es.num + gs*qs.num,
           iyr.flow = ry1*iyq.num + ry2*iyd.num + ry3*iyu.num + ry4*hy.num,
           isr.flow = rs1*isq.num + rs2*isd.num + rs3*isu.num + rs4*hs.num,
           sq.flow = qy*(i.num/s.num)*sy.num + qs*(i.num/s.num)*ss.num,
           se.flow = fy2*ynum + ((by1*iyd.num + by2*iyu.num)*ky + (bs1*isd.num + bs2*isu.num)*kys)*sy.num/ynum + fs2*snum +((bs1*isd.num + bs2*isu.num)*ks + (by1*iyd.num + by2*iyu.num)*kys)*ss.num/snum,
           eiy.flow = vy1*ey.num + gy*qy.num + vy2*ey.num,
           eis.flow = vs1*es.num + gs*qs.num + vs2*es.num,
           ei.flow = vy1*ey.num + gy*qy.num + vs1*es.num + gs*qs.num + vy2*ey.num + vs2*es.num,
           idq.flow = vy1*ey.num + gy*qy.num + vs1*es.num + gs*qs.num,
           eiu.flow = vy2*ey.num + vs2*es.num,
           r.flow = ry1*iyq.num + ry2*iyd.num + ry3*iyu.num + ry4*hy.num + rs1*isq.num + rs2*isd.num + rs3*isu.num + rs4*hs.num,
           h.flow = hy1*iyq.num +hy2*iyd.num + hs1*isq.num +hs2*isd.num,
           d.flow = wy*hy.num + ws*hs.num),
         i.num  = i.num,
         s.num = s.num,
         e.num = e.num,
         q.num = q.num,
         id.num = id.num,
         iu.num = iu.num,
         iq.num = iq.num,
         h.num = h.num,
         d.num = d.num,
         r.num = r.num,
         iy.incidence = sey.flow/sy.num,
         is.incidence = ses.flow/ss.num,
         i.incidence = (sey.flow + ses.flow)/(sy.num+ss.num),
         q.prev = q.num /(q.num+s.num),
         i.prev = i.num /num,
         h.prev = h.num/i.num,
         p.prev = 1-s.num/num)
  })
}


# sensitivity model (CV)
# a million sample size
# parameters should be in the order in the function;
# always reset the parameters from the default model

# imported case rate per days: 0.25: fy2 = 0.00000035  
# imported case rate per days: 1: fy2 = 0.0000014  
# imported case rate per days: 2: fy2 = 0.0000028  

# importing period: timpt = 5, timpt = 10, timpt = 30,

# contact old: ky = rep(10,5), ks = seq(3,15,3), kys = rep(3,5)
# contact old: ky = rep(15,5), ks = seq(3,15,3), kys = rep(3,5)
# contact young: ky = seq(3,15,3), ks = rep(7,5), kys = rep(3,5)
# contact both: ky = seq(3,15,3), ks = seq(3,15,3), kys = rep(3,5)

# self-quarantined: 0.5%: qy = 0.004, qs = 0.004
# self-quarantined: 5%: qy = 0.03, qs = 0.03
# self-quarantined: 10%: qy = 0.07, qs = 0.07

# asymptomatic cases: 10%: vy1 = 0.4, vs1 =0.55, vy2 = 0.05,vs2 = 0.075,  
# asymptomatic cases: 30%: vy1 = 0.22, vs1 =0.3, vy2 = 0.12,vs2 = 0.11,  
# asymptomatic cases: 60%: vy1 = 0.1, vs1 =0.18, vy2 = 0.2,vs2 = 0.3,  

# infectious of asym: 30% symptmatic: by2 = 0.014, bs2 = 0.03,
# infectious of asym: 70% symptmatic: by2 = 0.030, bs2 = 0.064,
# infectious of asym: 100% asymptmatic: by2 = 0.043, bs2 = 0.093,

# combined %asym and infectivity: 60% asym, 30% infectivity:
# vy1 = 0.1, vs1 =0.18, vy2 = 0.2,vs2 = 0.3, 
# by2 = 0.014, bs2 = 0.03,

param <- param.dcm(timpt = 20, 
                   fy1 = 0.00001, fs1 = 0, fy2 = 0.0000007,fs2 = 0, 
                   p = 0.07, qy = 0.01,qs = 0.01, 
                   uy = 0.00005, us = 0.00014,
                   ky = 10, ks = 7, kys = 3, 
                   by1 = 0.043, bs1= 0.093, by2 = 0.014, bs2 = 0.03, 
                   gy = 0.2, gs = 0.3, 
                   vy1 = 0.1, vs1 =0.18, vy2 = 0.2,vs2 = 0.3,
                   hy1 = 0.025, hs1 = 0.06, hy2 = 0.01, hs2 = 0.04, 
                   ry1 = 0.13, rs1 =0.07, ry2 = 0.13, rs2 = 0.07, 
                   ry3 = 0.18, rs3 = 0.10, ry4 = 0.133, rs4=0.047, 
                   wy = 0.007, ws = 0.012)



# 20% young people immuned, 20% elderly immuned
# Sensitivity analysis;

# size of region: 500,000: ynum = 400000,snum = 100000, sy.num = 320000, ss.num =80000 
# size of region: 100,000: ynum = 80000,snum = 20000, sy.num = 64000, ss.num =16000 
# size of region: 50,000: ynum = 40000,snum = 10000, sy.num = 32000, ss.num =8000 

# percent elderly: 10%: ynum = 900000,snum = 100000, sy.num = 720000, ss.num =80000
# percent elderly: 30%: ynum = 700000,snum = 300000, sy.num = 560000, ss.num =240000
# percent elderly: 40%: ynum = 600000,snum = 400000, sy.num = 480000, ss.num =320000

# % prior infection: 30%: ynum = 800000,snum = 200000, sy.num = 560000, ss.num =140000 
# % prior infection: 50%: ynum = 800000,snum = 200000, sy.num = 400000, ss.num =100000 
# % prior infection: 60%: ynum = 800000,snum = 200000, sy.num = 320000, ss.num =80000 
# % prior infection: 70%: ynum = 800000,snum = 200000, sy.num = 240000, ss.num =60000 

# order of parameters should be same as the diff equations, and same as the function output
# dNY,dNS,dQY,dQS,dSY,dSS,dEY,dES,dIYQ,dISQ,dIYD,dISD,dIYU,dISU,dHY,dHS, dRY,dRS,dDY,dDS,
init <- init.dcm(ynum = 800000,snum = 200000, 
                 qy.num = 0, qs.num=0, sy.num = 640000, ss.num =160000,
                 ey.num = 0,es.num=0, 
                 iyq.num = 0, isq.num = 0,iyd.num = 0,isd.num = 0, iyu.num = 0, isu.num = 0, 
                 hy.num = 0,hs.num = 0, ry.num = 0,rs.num = 0, dy.num = 0,ds.num = 0, 
                 ns.flow = 0, 
                 sqy.flow = 0,sqs.flow = 0, qys.flow = 0,qss.flow = 0, 
                 sey.flow = 0,ses.flow = 0,
                 ney.flow = 0,nes.flow = 0, 
                 eiyd.flow = 0,eisd.flow = 0, eiyu.flow = 0,eisu.flow = 0, 
                 qiyq.flow = 0,qisq.flow = 0,
                 iyqh.flow = 0,isqh.flow = 0, iyqr.flow = 0,isqr.flow = 0,
                 iydh.flow = 0,isdh.flow = 0, iydr.flow = 0,isdr.flow = 0,
                 iyur.flow = 0,isur.flow = 0,
                 hyr.flow = 0,hsr.flow = 0, hyd.flow = 0,hsd.flow = 0,
                 hy.flow = 0, hs.flow = 0, 
                 iy.flow = 0, is.flow = 0, 
                 iyr.flow = 0, isr.flow = 0,
                 sq.flow = 0, se.flow = 0, 
                 eiy.flow = 0, eis.flow = 0, ei.flow = 0, idq.flow = 0, eiu.flow = 0,
                 r.flow = 0, h.flow = 0, d.flow = 0)

control <- control.dcm(nsteps = 240, dt = 1, new.mod = SEIRQ)
cvmod <- dcm(param, init, control)

cvmod

cvmodepi<-as.data.frame(cvmod)

# plotting the model;

par(mfrow = c(1, 2))
plot(cvmod, y = "eiy.flow", main = "Age < 65",xlab ="Days",ylab = "Incident cases")
plot(cvmod, y = "eis.flow", main = "Age >= 65",xlab = "Days",ylab= "Incident cases",legend="full")

plot(cvmod, y = "iy.flow", main = "Age < 65",ylab="Incident symptomatic cases",xlab="Days")
plot(cvmod, y = "is.flow", main = "Age >= 65",ylab ="Incident symptomatic cases",xlab="Days",legend="full")

plot(cvmod, y = "eiyu.flow", main = "Age < 65",ylab="Incident asymptomatic cases",xlab="Days")
plot(cvmod, y = "eisu.flow", main = "Age >= 65",ylab ="Incident asymptomatic cases",xlab="Days",legend="full")

plot(cvmod, y = "hy.flow", main = "Age < 65",ylab="Hospitalizations",xlab="Days")
plot(cvmod, y = "hs.flow", main = "Age >= 65",ylab ="Hospitalizations",xlab="Days")

plot(cvmod, y = "hyd.flow", main = "Age < 65",ylab="Deaths",xlab="Days")
plot(cvmod, y = "hsd.flow", main = "Age >= 65",ylab ="Deaths",xlab="Days",legend="full")

# summarize results
peaktime<-cvmodepi$time[cvmodepi$ei.flow==max(cvmodepi$ei.flow)]
peaktime
cvmodepi$s.num[cvmodepi$time==peaktime]
cvmodepi$time[cvmodepi$ei.flow<= 10 & cvmodepi$time>100]

round(summary(cvmodepi$sq.flow))  # max all new quarantined
round(summary(cvmodepi$se.flow))  # max all new cases unquarantined
round(summary(cvmodepi$ei.flow))  # max all new cases
round(summary(cvmodepi$idq.flow))   # max symptomatic new cases
round(summary(cvmodepi$eiu.flow))   # max asymptomatic new cases
round(summary(cvmodepi$h.flow))   # max new hosps
round(summary(cvmodepi$d.flow))   # max new death

# max existing quaratined and susceptibles
max(cvmodepi$qy.num)
cvmodepi$sy.num[cvmodepi$qy.num==max(cvmodepi$qy.num)]
# max percent;
max(cvmodepi$qy.num/(cvmodepi$sy.num+cvmodepi$qy.num))

max(cvmodepi$qs.num)
cvmodepi$ss.num[cvmodepi$qs.num==max(cvmodepi$qs.num)]
max(cvmodepi$qs.num/(cvmodepi$ss.num+cvmodepi$qs.num))

length(cvmodepi$ei.flow)

cvmodepi$sq.sum <- cumsum(cvmodepi$sq.flow)
cvmodepi$se.sum <- cumsum(cvmodepi$se.flow)
cvmodepi$ei.sum <- cumsum(cvmodepi$ei.flow)
cvmodepi$idq.sum <- cumsum(cvmodepi$idq.flow)
cvmodepi$eiu.sum <- cumsum(cvmodepi$eiu.flow)
cvmodepi$h.sum <- cumsum(cvmodepi$h.flow)
cvmodepi$d.sum <- cumsum(cvmodepi$d.flow)

round(max(cvmodepi$se.sum))  # total all new cases unquarantined
# total cases, almost all population
round(c(max(cvmodepi$ei.sum),max(cvmodepi$idq.sum), max(cvmodepi$eiu.sum)))

round(c(max(cvmodepi$idq.flow),max(cvmodepi$eiu.flow),max(cvmodepi$h.flow),max(cvmodepi$d.flow),max(cvmodepi$h.sum),max(cvmodepi$d.sum)))


# young population
ypeaktime<-cvmodepi$time[cvmodepi$eiy.flow==max(cvmodepi$eiy.flow)]
ypeaktime
cvmodepi$sy.num[cvmodepi$time==ypeaktime]
cvmodepi$time[cvmodepi$eiy.flow<= 10 & cvmodepi$time>100]

round(summary(cvmodepi$sqy.flow))  # max all new quarantined
round(summary(cvmodepi$sey.flow))  # max all new cases unquaratined
round(summary(cvmodepi$eiy.flow))  # max all new cases
round(summary(cvmodepi$iy.flow))   # max symptomatic new cases
round(summary(cvmodepi$eiyu.flow))   # max asymptomatic new cases
round(summary(cvmodepi$hy.flow))   # max new hosps
round(summary(cvmodepi$hyd.flow))   # max new death

cvmodepi$sqy.sum<-cumsum(cvmodepi$sqy.flow)
cvmodepi$sey.sum<-cumsum(cvmodepi$sey.flow)
cvmodepi$eiy.sum<-cumsum(cvmodepi$eiy.flow)
cvmodepi$iy.sum<-cumsum(cvmodepi$iy.flow)
cvmodepi$eiyu.sum<-cumsum(cvmodepi$eiyu.flow)
cvmodepi$hy.sum<-cumsum(cvmodepi$hy.flow)
cvmodepi$hyd.sum<-cumsum(cvmodepi$hyd.flow)

round(max(cvmodepi$sey.sum))  # max all new cases from unquarantined

# total cases, almost all population
round(c(max(cvmodepi$eiy.sum),max(cvmodepi$iy.sum), max(cvmodepi$eiyu.sum)))

round(c(max(cvmodepi$iy.flow),max(cvmodepi$eiyu.flow),max(cvmodepi$hy.flow),max(cvmodepi$hyd.flow),max(cvmodepi$hy.sum),max(cvmodepi$hyd.sum)))


#old population
speaktime<-cvmodepi$time[cvmodepi$eis.flow==max(cvmodepi$eis.flow)]
speaktime
cvmodepi$ss.num[cvmodepi$time==speaktime]

round(summary(cvmodepi$sqs.flow))  # max all new quarantined
round(summary(cvmodepi$ses.flow))  # max all new cases from unquarantined
round(summary(cvmodepi$eis.flow))   # max  new cases
round(summary(cvmodepi$is.flow))   # max symptomatic new cases
round(summary(cvmodepi$eisu.flow))   # max asymptomatic new cases
round(summary(cvmodepi$hs.flow))   # max new hosps
round(summary(cvmodepi$hsd.flow))   # max new death

cvmodepi$sqs.sum<-cumsum(cvmodepi$sqs.flow)
cvmodepi$ses.sum<-cumsum(cvmodepi$ses.flow)
cvmodepi$eis.sum<-cumsum(cvmodepi$eis.flow)
cvmodepi$is.sum<-cumsum(cvmodepi$is.flow)
cvmodepi$eisu.sum<-cumsum(cvmodepi$eisu.flow)
cvmodepi$hs.sum<-cumsum(cvmodepi$hs.flow)
cvmodepi$hsd.sum<-cumsum(cvmodepi$hsd.flow)

round(max(cvmodepi$ses.sum))  # total all new cases from unquanratined

# total cases, almost all population
round(c(max(cvmodepi$eis.sum),max(cvmodepi$is.sum), max(cvmodepi$eisu.sum)))

round(c(max(cvmodepi$is.flow),max(cvmodepi$eisu.flow),max(cvmodepi$hs.flow),max(cvmodepi$hsd.flow),max(cvmodepi$hs.sum),max(cvmodepi$hsd.sum)))

# final results:
# total
peaktime
cvmodepi$time[cvmodepi$ei.flow<= 5 & cvmodepi$time>100]
round(c(max(cvmodepi$idq.flow),max(cvmodepi$eiu.flow),max(cvmodepi$h.flow),max(cvmodepi$d.flow),max(cvmodepi$h.sum),max(cvmodepi$d.sum)))

# young
ypeaktime
cvmodepi$time[cvmodepi$eiy.flow<= 5 & cvmodepi$time>100]
round(c(max(cvmodepi$iy.flow),max(cvmodepi$eiyu.flow),max(cvmodepi$hy.flow),max(cvmodepi$hyd.flow),max(cvmodepi$hy.sum),max(cvmodepi$hyd.sum)))

# old
speaktime
cvmodepi$time[cvmodepi$eis.flow<= 5 & cvmodepi$time>100]
round(c(max(cvmodepi$is.flow),max(cvmodepi$eisu.flow),max(cvmodepi$hs.flow),max(cvmodepi$hsd.flow),max(cvmodepi$hs.sum),max(cvmodepi$hsd.sum)))










# write out analysis for asymptomatic cases

asymmod <- data.frame(defhy=defmodepi$hy.flow,
                      defhs=defmodepi$hs.flow,
                      defh =defmodepi$h.flow,
                      defeiy =defmodepi$eiy.flow,
                      defeis =defmodepi$eis.flow,
                      defei =defmodepi$ei.flow,
                      defeiyu =defmodepi$eiyu.flow,
                      asymhy=cvmodepi$hy.flow,
                      asymhs=cvmodepi$hs.flow,
                      asymh = cvmodepi$h.flow,
                      asymeiy =cvmodepi$eiy.flow,
                      asymeis =cvmodepi$eis.flow,
                      asymei =cvmodepi$ei.flow,
                      asymeiyu =cvmodepi$eiyu.flow,
                      time=defmodepi$time)

write.dta(asymmod,"asymmod2.dta")
