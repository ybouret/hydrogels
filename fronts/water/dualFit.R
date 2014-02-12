pix2pos = 5.0/427.0;
pix2tmx = 10.0;

#first front times
t1 = c( 
 3,  
 8 , 
 14 ,
 20 ,
 27 ,
 32 ,
 43 ,
 55 ,
 64 ,
 72 ,
 83 ,
 89 ) * pix2tmx;
 
 #first front position
 x1 = c(
 377, 
 428, 
 481, 
 529, 
 574, 
 626, 
 684, 
 746, 
 794, 
 829, 
 883, 
 929 
 ) * pix2pos;
 
z1 = x1^2;

#plot(t1,z1,col='blue');
fit1 = lm(z1~t1);
summary(fit1)
b1=fit1$coeff[1];
a1=fit1$coeff[2];
#lines(t1,a1*t1+b1,col='blue');

#second front times
t2 = c(
6,   
12,  
20 , 
27  ,
36  ,
46  ,
58  ,
67  ,
79  ,
92  ,
104 ) * pix2tmx;

# second positions
x2 = c(
444,   
500,   
569,   
629,   
689,   
750,   
818,   
869,   
929,   
995,  
1043  
)* pix2pos;

z2 = x2^2;

#points(t2,z2,col='red');
fit2 = lm(z2~t2);
summary(fit2)

b2=fit2$coeff[1];
a2=fit2$coeff[2];
#lines(t2,a2*t2+b2,col='red');


#joining
z_all=c(z1,z2);

n1 = length(t1);
n2 = length(t2);

t1_all=c(t1,rep(0,n2));

t2_all=c(rep(0,n1),t2);

fit_all=lm(z_all~t1_all+t2_all);

z0=fit_all$coeff[1];
slope1=fit_all$coeff[2];
slope2=fit_all$coeff[3];

plot(t1,x1,col='blue',xlim=c(min(c(t1,t2)),max(c(t1,t2))), ylim=c(min(c(x1,x2)),max(c(x1,x2))));
points(t2,x2,col='red');
lines(t1,sqrt(z0+slope1*t1),col='blue');
lines(t2,sqrt(z0+slope2*t2),col='red');


 