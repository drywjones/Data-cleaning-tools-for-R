## Outlier correcting function with absolute or relative metrics
## for cutoffs (extreme differnces from previous values) values 
## in time series data with periodic patterns
data.cleaner<-function(data.vect,
                             crit.num,
                             trails=24,
                             rel.or.abs="absolute",
                             differential=1,
                             center.dat=0.95,
                             TIMESTAMP,
                             flex=9,
                             data.type="periodic"){
  require(splines)
  require(lubridate)
  ## data.vect    - these are the observations that including outliers. This should be in 
  ##                vector form.
  ## crit.num     - this is the critical value above which data is deemed to be an 
  ##                outlier (for absolute cutoff), or the number of standard deviations
  ##                from the central 95% of data for relative cutoffs, the smaller the central
  ##                data range (eg. 90%, 85% etc) the smaller this value should be.
  ## trails       - how many of the trailing observations should be included in the 
  ##                spline fit used to replace outlier? Should use as many observations
  ##                as make up at least one cycle (if periodic) otherwise non-periodic will
  ##                be the assumption for spline fits. 
  ## rel.or.abs   - switch that determines whether a relative cutoff (number of standard
  ##                deviations), or an absolute cutoff (critical number) is used.
  ## differential - determines how much more (or less) sensitive the outliers are to 
  ##                negative values.
  ## center.dat   - determines the central percent of data that is used for relative cutoff 
  ##                estimates. 
  ## TIMESTAMP    - these are the times of observations - if no time vector is provided
  ##                then a sequence periodic values based on trails (if periodic) or a 
  ##                consecutive sequence of numbers based on total length of the data
  ##                vector will be created.
  if(missing(TIMESTAMP)==TRUE){
    ## create periodic data for spline fitting based on trails
    if(data.type=="periodic"){
      time=rep(seq(1:trails),ceiling(length(data.vec)/trails))[1:length(data.vect)]
    }else{
      time=seq(1:length(data.vect))
    }
    
    
  }else{time=TIMESTAMP}
    if(rel.or.abs=="absolute"){
      
      crit.val.pos<-crit.num
      }else{
        cent.dif<-1-center.dat
        bot.q<-cent.dif/2
        top.q<-1-cent.dif/2
        quants<-quantile(data.vect,c(bot.q,top.q),na.rm=T)
        data.vect.cent<-data.vect[data.vect>=quants[1]&data.vect<=quants[2]]
        crit.val.pos<-crit.num*sd(data.vect.cent,na.rm=T)
      }
  
  crit.val.neg<--1*crit.val.pos/differential
  num.iter<-length(data.vect)
  orig.vect<-data.vect
  for(i in 1:num.iter){
    if(i==1){
      orig.vect[i]<-orig.vect[i]
    }else{
      if(is.na(orig.vect[(i-1)])&&is.na(orig.vect[(i+1)])){
        
        orig.vect[i]<-orig.vect[i]
      }else{
        if(is.na(orig.vect[i])){
          grow.amt<-mean(orig.vect[(i-1):(i+1)],na.rm=T)
        }else{
          grow.amt<-orig.vect[i]-mean(c(orig.vect[(i-1)],orig.vect[(i+1)]),na.rm=T)
        }
       
        if(grow.amt>=crit.val.pos|grow.amt<=crit.val.neg){
          if(i<=1.5*flex){
            orig.vect[i]<-median(orig.vect[1:(i+1)],na.rm=T)
          }
          else{
          if(i<=(trails+1)){
            orig.vect[i]<-median(orig.vect[(i-3):(i+3)],na.rm=T)
          }else{
            if(length(orig.vect[!is.na(orig.vect[c((i-trails):(i-1),i-1)])])<=ceiling(trails/flex+1)){
              orig.vect[i]<-orig.vect[i]
            }else{
            
            a<-time[c((i-trails):i)]
            b<-c(orig.vect[c((i-trails):(i-1),i-1)])
            dat<-data.frame(a,b)
            preds<-(predict(lm(b~bs(a,df=(ceiling(trails/flex+1))),data=dat,na.action=na.omit)))
            preds<-preds[!is.na(preds)]
            if(length(preds)==0){
              orig.vect[i]<-orig.vect[i]
            }else{
              orig.vect[i]<-preds[length(preds)]
            }
            
          }
         }
        }
       }
      }
    }
  }
  clean.dat<-orig.vect
  return(clean.dat)}
#
# Create data that somewhat represents general patterns over time 
# based on observed air temperature data at a location on the planet 
# in 2010.

# create hourly time series data:
times<-seq(ISOdate(2010,1,1),ISOdate(2011,1,1),"hour")
# create an hour of the day vector from time series:
hours<-as.numeric(format(times,format="%H"))
# create a calendar day vector from the time series:
cal.day<-yday(times)
# create data from modeled trends with added random noise:
orig.data<-sample(1:100,length(hours),replace=T)/100+sin((hours-9)*pi/12)+(4+3.528e-02*cal.day-9.882e-04*cal.day^2+1.523e-05*cal.day^3-6.562e-08*cal.day^4+ 8.367e-11*cal.day^5)
# Duplicate this data to create messy data
data<-orig.data
# Insert some missing data points at random:
data[c(40:100,sample(1:length(hours),500))]<-NA
# Create indexes for outliers that will be too high
high.sample<-sample(1:length(hours),2500)
# Create indexes for outliers that will be too low
low.sample<-sample(1:length(hours),3500)
# Add noise to data for high and low outliers:
data[low.sample]<-data[low.sample]-6*sample(1:100,length(low.sample),replace=T)/100
data[high.sample]<-data[high.sample]+8*sample(1:100,length(high.sample),replace=T)/100
# Make a data.frame from the derived data:
data.fr<-data.frame(times=times,hours=hours,data=data,orig.data=orig.data)

## Now create cleaner data - using hours for time series in this case as we 
## are assuming that daily temperature patterns are similar from one day to 
## another and that periodicity is reflected by the hour of day rather than 
## calendar date and time.
data.fr$clean.data<-data.cleaner(data.fr$data,
                                       crit.num=1,
                                       trails=24,
                                       rel.or.abs="absolute",
                                       differential=1,
                                       TIMESTAMP=data.fr$hours,
                                       flex=12)

## same thing but adjust for differences in negative versus positive drops:
  data.fr$clean.data.low.hi<-data.cleaner(data.fr$data,
                                         crit.num=1,
                                         trails=24,
                                         "absolute",
                                         differential=1.33,
                                         TIMESTAMP=data.fr$hours,
                                         flex=12)
  ## using a relative cutoff instead
  data.fr$clean.data.rel<-data.cleaner(data.fr$data,
                                         crit.num=.25,
                                         trails=24,
                                         rel.or.abs="relative",
                                         differential=1.33,
                                         TIMESTAMP=data.fr$hours,
                                         center.dat=.95,
                                         flex=12)
## Check on differences in means, quartiles, and extremes:  
  summary(data.fr$orig.data-data.fr$data)
  summary(data.fr$orig.data-data.fr$clean.data)
  summary(data.fr$orig.data-data.fr$clean.data.low.hi)
  summary(data.fr$orig.data-data.fr$clean.data.rel)
  
## Visualizing the impact on data points:  
library(ggplot2)

## Create a plot fo the whole year:
year.plot<-(ggplot(data.fr,aes(times,data),color="black")+
              geom_point(aes(times,data),color="red",alpha=.25)+
              geom_point(aes(times,clean.data),color="blue",alpha=.25)+
              geom_point(aes(times,clean.data.rel),color="grey",alpha=.25)+
              geom_point(aes(times,orig.data),color="green",alpha=.25)
)
print(year.plot)

## zoomed in view to small subset of data:
sub.year<-data.fr[1300:2300,]
sub.year.plot<-(ggplot(sub.year,aes(times,data))+
                             geom_point(aes(times,data),color="red")+
                             geom_point(aes(times,clean.data),color="blue")+
                             geom_point(aes(times,orig.data),color="grey")+
                             geom_point(aes(times,clean.data.rel),color="green")
)
print(sub.year.plot)
  


