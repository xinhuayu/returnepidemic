library(EpiModel)
# program for contacts

# customized SEIR with quarantine and imported cases;

SEIRQ <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # aggregate parameters over young and old
    num = ynum + snum
    q.num = qy.num + qs.num
    s.num = sy.num + ss.num
    e.num = ey.num + es.num
    #total symptomatic prev cases;
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

# contact old: ky = rep(15,5), ks = seq(3,15,3), kys = rep(3,5)
# contact old: ky = rep(10,5), ks = seq(3,15,3), kys = rep(3,5)
# contact old: ky = rep(7,5), ks = seq(3,15,3), kys = rep(3,5)
# contact young: ky = seq(3,15,3), ks = rep(7,5), kys = rep(3,5)
# contact both: ky = seq(3,15,3), ks = seq(3,15,3), kys = rep(3,5)

# extreme case: 
# high risk:
# 1 case per day (fy2 = 0.0000014) continuous (timp = 200 days)
# 30% asymptomatic and same infectivity of asymptomatic: 
# vy1 = 0.18, vs1 =0.22, vy2 = 0.13,vs2 = 0.15,
# by2 = 0.043, bs2= 0.093, 


# Low risk: Ideal situtation
# 1 case per two days (fy2 = 0.0000007) continuous (timp = 20 days)
# 60% percent of asymptomatic, 30% infectivity 
# vy1 = 0.1, vs1 =0.18, vy2 = 0.2,vs2 = 0.3, 
# by2 = 0.014, bs2 = 0.03,

param <- param.dcm(timpt = 200, 
                   fy1 = 0.00001, fs1 = 0, fy2 = 0.0000014,fs2 = 0, 
                   p = 0.07, qy = 0.01,qs = 0.01, 
                   uy = 0.00005, us = 0.00014,
                   ky = rep(10,5), ks = seq(3,15,3), kys = rep(3,5), 
                   by1 = 0.043, bs1= 0.093, by2 = 0.043, bs2= 0.093, 
                   gy = 0.2, gs = 0.3, 
                   vy1 = 0.18, vs1 =0.22, vy2 = 0.13,vs2 = 0.15, 
                   hy1 = 0.025, hs1 = 0.06, hy2 = 0.01, hs2 = 0.04, 
                   ry1 = 0.13, rs1 =0.07, ry2 = 0.13, rs2 = 0.07, 
                   ry3 = 0.18, rs3 = 0.10, ry4 = 0.133, rs4=0.047, 
                   wy = 0.007, ws = 0.012)

# 20% young people immuned, 20% elderly immuned
pyzero = 0.2
pszero = 0.2

# dNY,dNS,dSY,dSS,dQY,dQS,dEY,dES,dIYQ,dISQ,dIYD,dISD,dIYU,dISU,dHY,dHS, dRY,dRS,dDY,dDS,


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
# changing names depending on the paramters;
par(mfrow = c(1, 2))

# names(cvmod$epi$eiy.flow)<-c("3 contacts","6 contacts","9 contacts","12 contacts","15 contacts")
# names(cvmod$epi$eiy.flow)<-c("15 contacts","15 contacts","15 contacts","15 contacts","15 contacts")
 names(cvmod$epi$eiy.flow)<-c("10 contacts","10 contacts","10 contacts","10 contacts","10 contacts")

 names(cvmod$epi$eis.flow)<-c("3 contacts","6 contacts","9 contacts","12 contacts","15 contacts")
# names(cvmod$epi$eis.flow)<-c("7 contacts","7 contacts","7 contacts","7 contacts","7 contacts")

plot(cvmod, y = "eiy.flow", main = "Age < 65",xlab ="Days",ylab = "Incident cases",legend="full")
plot(cvmod, y = "eis.flow", main = "Age >= 65",xlab = "Days",ylab= "Incident cases",legend="full")


# names(cvmod$epi$iy.flow)<-c("3 contacts","6 contacts","9 contacts","12 contacts","15 contacts")
# names(cvmod$epi$iy.flow)<-c("15 contacts","15 contacts","15 contacts","15 contacts","15 contacts")
 names(cvmod$epi$iy.flow)<-c("10 contacts","10 contacts","10 contacts","10 contacts","10 contacts")

 names(cvmod$epi$is.flow)<-c("3 contacts","6 contacts","9 contacts","12 contacts","15 contacts")
# names(cvmod$epi$is.flow)<-c("7 contacts","7 contacts","7 contacts","7 contacts","7 contacts")

plot(cvmod, y = "iy.flow", main = "Age < 65",ylab="Incident symptomatic cases",xlab="Days",legend="full")
plot(cvmod, y = "is.flow", main = "Age >= 65",ylab ="Incident symptomatic cases",xlab="Days",legend="full")

# names(cvmod$epi$eiyu.flow)<-c("3 contacts","6 contacts","9 contacts","12 contacts","15 contacts")
# names(cvmod$epi$eiyu.flow)<-c("15 contacts","15 contacts","15 contacts","15 contacts","15 contacts")
 names(cvmod$epi$eiyu.flow)<-c("10 contacts","10 contacts","10 contacts","10 contacts","10 contacts")

 names(cvmod$epi$eisu.flow)<-c("3 contacts","6 contacts","9 contacts","12 contacts","15 contacts")
# names(cvmod$epi$eisu.flow)<-c("7 contacts","7 contacts","7 contacts","7 contacts","7 contacts")

plot(cvmod, y = "eiyu.flow", main = "Age < 65",ylab="Incident asymptomatic cases",xlab="Days",legend="full")
plot(cvmod, y = "eisu.flow", main = "Age >= 65",ylab ="Incident asymptomatic cases",xlab="Days",legend="full")

# names(cvmod$epi$hy.flow)<-c("3 contacts","6 contacts","9 contacts","12 contacts","15 contacts")
# names(cvmod$epi$hy.flow)<-c("15 contacts","15 contacts","15 contacts","15 contacts","15 contacts")
 names(cvmod$epi$hy.flow)<-c("10 contacts","10 contacts","10 contacts","10 contacts","10 contacts")

 names(cvmod$epi$hs.flow)<-c("3 contacts","6 contacts","9 contacts","12 contacts","15 contacts")
# names(cvmod$epi$hs.flow)<-c("7 contacts","7 contacts","7 contacts","7 contacts","7 contacts")

plot(cvmod, y = "hy.flow", main = "Age < 65",ylab="Hospitalizations",xlab="Days",legend="full")
plot(cvmod, y = "hs.flow", main = "Age >= 65",ylab ="Hospitalizations",xlab="Days",legend="full")

# names(cvmod$epi$hyd.flow)<-c("3 contacts","6 contacts","9 contacts","12 contacts","15 contacts")
# names(cvmod$epi$hyd.flow)<-c("15 contacts","15 contacts","15 contacts","15 contacts","15 contacts")
 names(cvmod$epi$hyd.flow)<-c("10 contacts","10 contacts","10 contacts","10 contacts","10 contacts")

 names(cvmod$epi$hsd.flow)<-c("3 contacts","6 contacts","9 contacts","12 contacts","15 contacts")
# names(cvmod$epi$hsd.flow)<-c("7 contacts","7 contacts","7 contacts","7 contacts","7 contacts")

plot(cvmod, y = "hyd.flow", main = "Age < 65",ylab="Deaths",xlab="Days",legend="full")
plot(cvmod, y = "hsd.flow", main = "Age >= 65",ylab ="Deaths",xlab="Days",legend="full")


# summarize results
# peak time is visual
# max hosp

round(aggregate(cvmodepi$hy.flow,by=list(cvmodepi$run), summary))   # max new hosps
round(aggregate(cvmodepi$hyd.flow,by=list(cvmodepi$run), summary))   # max new deaths
round(aggregate(cvmodepi$eiy.flow,by=list(cvmodepi$run), summary))   # max new cases
round(aggregate(cvmodepi$iy.flow,by=list(cvmodepi$run), summary))   # max new symptomatic cases
round(aggregate(cvmodepi$eiyu.flow,by=list(cvmodepi$run), summary))   # max new asymptomatic cases

round(aggregate(cvmodepi$hs.flow,by=list(cvmodepi$run), summary))   # max new hosps
round(aggregate(cvmodepi$hsd.flow,by=list(cvmodepi$run), summary))   # max new deaths
round(aggregate(cvmodepi$eis.flow,by=list(cvmodepi$run), summary))   # max new cases
round(aggregate(cvmodepi$is.flow,by=list(cvmodepi$run), summary))   # max new symptomatic cases
round(aggregate(cvmodepi$eisu.flow,by=list(cvmodepi$run), summary))   # max new asymptomatic cases


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
