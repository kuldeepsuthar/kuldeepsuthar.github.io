function maincode(P0,d,b,tf,tw,fy,axs,t1,t4,t5){
	  var r =0;
	  var t2 = t1;
	  var t3 = t1;
	  var t6 = t5;
	  var t7 = t5;
	  var gs1 = 5;
	  var gs2 = 2;
	  var gs3 = 20;
	  var gs4 = 10;
	  var S = [];
	  for (j=0; j<91; j++){
		    S[j] = 0;
	  }
	  var cur = [];
	  [matrix,x,y, Gx, Gy,dx,dy,Area,area,Ix,Iy,ix,iy,resi] = Sectionproperties(b,d,tf,tw,gs1,gs2,gs3,gs4,r);
	  [temp] = temperature(matrix,b,d,tf,gs1,gs2,gs3,dx,dy,Gx,Gy,x,y,t1,t2,t3,t4,t5,t6,t7);
	  [tstrain] = tempelong(temp,x,y,matrix);
	  var load = P0;
	  P0 = load*1000;
	  console.log(P0);
	  var Py = area*fy;
	  
	  console.log(Py);
  var E = 200000;
    var ststrain=0;
  for (var i=0; i<=x; i++){
    for (var j = 0; j <=y; j++){
      if (matrix[i][j]==1){
        ststrain = ststrain+ ((tstrain[i][j]*Area[i][j])/area);
      }
    }
  }
  var p = (P0/(area*E))-ststrain;
	  var phi = 0;
	  var k = -1;
	  for (ki=-45; ki<46; ki++){
		
		k=k+1;
		
		if (ki<-30){
			phi = (-0.00003+(0.00000525*5*(ki+29)));
		}
		if (ki<-20 && ki>=-30){
			phi = (-0.00003+(0.00000525*(ki+21)));
		}
		if (ki<-10 && ki>=-20){
			phi = (-0.0000125+(0.00000175*(ki+10)));
		}
		if (ki<0 && ki>=-10){
			phi = (0.00000125*ki);                      
		}
		if (ki>=0 && ki<=10){
			phi = (0.00000125*ki);                      
		}
		if (ki>10 && ki<=20){
			phi = (0.0000125+(0.00000175*(ki-10)));       
		}
		if (ki>20 && ki<=30){
			phi = (0.00003+(0.00000525*(ki-21)));       
		}
		if (ki>30){
			phi = (0.00003+(0.00000525*5*(ki-29)));         
		}
		phi = phi*500;
		[Mx] = Momentcal(dx,dy,Area,matrix,Ix,Iy,phi,x,y,axs,P0,area,b,d,Gx,Gy,resi,tstrain,temp,fy,p);
	  var q = Mx.length;
	  cur[k] = phi;
	  for (j=0; j<q; j++){
	  S[k] = S[k]+Mx[j];
	  }
	  }
	  return [cur,S];
	  
	}
	function maincodeminor(P0,d,b,tf,tw,fy,axs,t1,t2,t3){
	  var r =0;
	  var t5 = t1;
	  var t4 = t2;
	  var t6 = t2;
	  var t7 = t3;
	  var gs1 = 4;
	  var gs2 = 4;
	  var gs3 = 15;
	  var gs4 = 15;
	  var S = [];
	  for (j=0; j<91; j++){
		S[j] = 0;
	  }
	  var cur = [];
	  [matrix,x,y, Gx, Gy,dx,dy,Area,area,Ix,Iy,ix,iy,resi] = Sectionproperties(b,d,tf,tw,gs1,gs2,gs3,gs4,r);
	  [temp] = temperature(matrix,b,d,tf,gs1,gs2,gs3,dx,dy,Gx,Gy,x,y,t1,t2,t3,t4,t5,t6,t7);
	  [tstrain] = tempelong(temp,x,y,matrix);
        var ststrain=0;
      var E = 200000;
  for (var i=0; i<=x; i++){
    for (var j = 0; j <=y; j++){
      if (matrix[i][j]==1){
        ststrain = ststrain+ (tstrain[i][j]*Area[i][j]/area);
      }
    }
  }
  var p = (P0/(area*E))-ststrain;
  var load = P0;
	  P0 = load*1000;
      console.log(P0);
	  var Py = area*fy;
	  console.log(Py);
	  var phi = 0;
	  var k = -1;
	  for (ki=-45; ki<46; ki++){
		
		k=k+1;
		
		if (ki<-30){
			phi = (-0.00003+(0.00000525*5*(ki+29)));
		}
		if (ki<-20 && ki>=-30){
			phi = (-0.00003+(0.00000525*(ki+21)));
		}
		if (ki<-10 && ki>=-20){
			phi = (-0.0000125+(0.00000175*(ki+10)));
		}
		if (ki<0 && ki>=-10){
			phi = (0.00000125*ki);                      
		}
		if (ki>=0 && ki<=10){
			phi = (0.00000125*ki);                      
		}
		if (ki>10 && ki<=20){
			phi = (0.0000125+(0.00000175*(ki-10)));       
		}
		if (ki>20 && ki<=30){
			phi = (0.00003+(0.00000525*(ki-21)));       
		}
		if (ki>30){
			phi = (0.00003+(0.00000525*5*(ki-29)));         
		}
		phi = phi*500;
		[Mx] = Momentcal(dx,dy,Area,matrix,Ix,Iy,phi,x,y,axs,P0,area,b,d,Gx,Gy,resi,tstrain,temp,fy,p);
	  var q = Mx.length;
	  cur[k] = phi;
	  for (j=0; j<q; j++){
	  S[k] = S[k]+Mx[j];
	  }
	  }
	  return [cur,S];
	}

function Momentcal(dx,dy,Area,matrix,Ix,Iy,phi,x,y,axs,P0,area,b,d,Gx,Gy,resi,tstrain,temp,fy,p){
  var E = 200000;
  [strainCG] = Newtonrap(phi,dx,dy,Gx,Gy,P0,x,y,axs,Area,matrix,p,b,d,resi,tstrain,temp,E,fy,area);
  var Yna = Gy-(strainCG/phi);
  var Xna = Gx-(strainCG/phi);
  
  [Mx] = momentstep(x,y,matrix,dx,dy,Area,Ix,Iy,Gx,Gy,axs,phi,strainCG,b,d,resi,tstrain,temp,E,fy);
  return [Mx];
}
function momentstep(x,y,matrix,dx,dy,Area,Ix,Iy,Gx,Gy,axs,phi,strainCG,b,d,resi,tstrain,temp,E,fy){
var Mx = [];
var v = 0;
var M = 0;
var strain;
for(var i=0; i<=x; i++) {
    Mx[i] = [];
}
if (axs==0){
    for (var i=0; i<=x; i++){
        v =0;
        for (var j=0; j<=y; j++){
            if (matrix[i][j] ==1){
                strain = strainCG+((phi*(dy[i][j]-(d-Gy)))/1000)+resi[i][j]+tstrain[i][j];
              
              [f,slope]=stress(strain,i,j,temp,E,fy);
                M = ((f*Area[i][j]*((dy[i][j]-(d-Gy))))/1000);
                v= v+M;
            }
        }
        Mx[i] = v;
    }
}
if (axs==1){
    for (var i=0; i<=x; i++){
        var v =0;
        for (var j=0; j<=y; j++){
            if (matrix[i][j] ==1){
                
                strain = strainCG+((phi*(dx[i][j]-(b-Gx)))/1000)+resi[i][j]+tstrain[i][j];
                
              [f,slope]=stress(strain,i,j,temp,E,fy);
                M = ((f*Area[i][j]*((dx[i][j]-(b-Gx))))/1000);
                v= v+M;
              
            }
         
        }
      
        Mx[i] = v;
    }
}
return [Mx];
  
}
function Newtonrap(phi,dx,dy,Gx,Gy,P0,x,y,axs,Area,matrix,p,b,d,resi,tstrain,temp,E,fy,area){
var error = 10000;
var countr = 0;
var z = p;

if (P0!=0){
    var tolerr = Math.abs(0.001*P0);
}
else{
    var tolerr= 0.1;
}
while (Math.abs(error)>tolerr){
    
    var P = 0;
    var ww = 0;
    var www = 0;
    var forcevar = 0;
    var strainCG = z;
    var strain =0;
    var k = 0;
    countr = countr +1;
    for (var i = 0; i <= x; i++){
        for (var j = 0; j <= x; j++){
            if (matrix[i][j] ==1){
                if (axs ==0){
                   
                    strain = strainCG+((phi*(dy[i][j]-(d-Gy)))/1000)+resi[i][j]+tstrain[i][j];
                }
                if (axs ==1){
                    
                    strain = strainCG+((phi*(dx[i][j]-(b-Gx)))/1000)+resi[i][j]+tstrain[i][j];
                }
                
                [f,slope]=stress(strain,i,j,temp,E,fy);
                k = f*Area[i][j];
                www = P;
                P = k+www;
                
                //ww = forcevar;
                // forcevar = ww+((slope*Area[i][j]));
            }
        }
    }
 // if (countr<5){
 //         forcevar= area*200000*0.01; 
//  }
  if (countr > 0 && countr < 15){
    forcevar = area*200000*0.07;
  }
  if (countr > 14 && countr < 30){
    forcevar = area*2000000*0.1;
  }
    if (countr > 29 && countr < 40){
    forcevar = area*200000*0.5;
  }
      if (countr > 39 && countr < 80){
    forcevar = area*200000;
  }
      if (countr > 79 && countr < 100){
    forcevar = area*200000*5;
  }
  if (countr > 99 && countr < 120){
    forcevar = area*200000*7;
  }
  if (countr > 119 && countr < 150){
    forcevar = area*200000*25;
  }
    if (countr > 149){
    forcevar = area*200000*50;
  }
  


  
    error = P0-P;
         if (error>0){
             z = strainCG+((error/forcevar));
            
         }
         if (error<0){
             z = strainCG+((error/forcevar));
            
         }
    
         error = Math.abs(error);
     if (countr ==200){
       var per_error = error*100/P0;
       console.log(per_error);
       break;
     }
}
strainCG = z;
return [strainCG];
}

function Sectionproperties(b,d,tf,tw,gs1,gs2,gs3,gs4,r){
var x1 = (tf/gs1);                   
var x2 = (tw/gs2);                  
var x3 = (d-(2*tf))/gs3;             
var x4 = ((b-tw)/(2*gs4));          
var x = (2*gs1)+gs3-1;                
var y = gs2+(2*gs4)-1;                
var Gx = (b/2);                  
var Gy = (d/2);
  var resi = [];
  var matrix = [];
  var dx = [];
  var dy = [];
  var Area = [];
  var Ix = [];
  var Iy = [];
  var area = 0;
for(var i=0; i<=x; i++) {
    resi[i]=[];
    matrix[i] = [];
    dx[i] = [];
    dy[i] = [];
    Area[i] = [];
    Ix[i] = [];
    Iy[i] = [];
    for(var j=0; j<=y; j++) {
      resi[i][j] = 0;
      matrix[i][j] = 0;
      dx[i][j] = 0;
      dy[i][j] = 0;
      Area[i][j] = 0;
      Ix[i][j] = 0;
      Iy[i][j] = 0;
      
    }
}
  for (var i = 0; i <= x; i++){
    for (var j = 0; j <= y; j++){
             if (i>(gs1-1) && i<=(x-gs1) && (j <=gs4-1 ||   j >gs2+gs4-1 )){
                
                matrix[i][j] = 0;
               
             }
             else{
		        matrix[i][j] = 1;
	         }
    }
  }
for (var i = 0; i <= x; i++){
    for (var j = 0; j <= y; j++){
      if (matrix[i][j]==1){
             if (j<=gs4-1){
                dx[i][j] = Math.abs((((j+0.5)*x4)));
             }
             if (j>(y-gs4)){
                dx[i][j]= Math.abs(b-((y-j+0.5)*x4));
             }
             if (j>gs4-1 && j<=(gs2+gs4)-1){
                dx[i][j] = Math.abs(((((j-gs4+0.5)*x2))+((b-tw)/2)));

             }
             if (i<=gs1-1){
                dy[i][j] = Math.abs((d)-((i+0.5)*x1));
             }
             if (i>(gs3+gs1)-1){
                dy[i][j] = Math.abs((tf)-((i-gs1-gs3+0.5)*x1));
             }
             if (i>(gs1-1) && i<=(x-gs1)){
                dy[i][j] = Math.abs((d)-((((i-gs1+0.5)*x3))+tf));
             }
      }
    }

}
  for (var i = 0; i <= x; i++){
    for (var j = 0; j <= y; j++){
      if (matrix[i][j]==1){
        if ((i<=gs1-1||(i>(gs1+gs3)-1 && i<=x)) && (j<=gs4-1 || ( j> (gs2+gs4)-1 && j<=y))){
           if (dx[i][j]< b/2){
             resi[i][j] = r-(2*r*dx[i][j]*2/b);
           }else{
             resi[i][j] = -r+(2*r*(dx[i][j]-0.5*b)*2/b);
           }

          Area[i][j] = x1*x4;
          Ix[i][j] = (x4*(x1**3))/12;
          Iy[i][j] = (x1*(x4**3))/12;
        }
        if ((i<=gs1-1 ||(i>(gs1+gs3)-1 && i<=x) ) && j>gs4-1 && j<=(gs2+gs4)-1){
          if (dx[i][j]< b/2){
             resi[i][j] = r-(2*r*dx[i][j]*2/b);
          }else{
                 resi[i][j] = -r+(2*r*(dx[i][j]-0.5*b)*2/b);
          }
        
          Area[i][j] = x1*x2;
          Ix[i][j] = (x2*(x1**3))/12;
          Iy[i][j] = (x1*(x2**3))/12;
        }
        if (i>gs1-1 && i<=(gs1+gs3)-1 && (j<=gs4-1 ||(j> (gs2+gs4)-1 && j<=y))){
          Area[i][j] = x3*x4;
          Ix[i][j] = (x4*(x3**3))/12;
          Iy[i][j] = (x3*(x4**3))/12;
        }
        if (i>gs1-1 && i<=(gs1+gs3)-1 && j>gs4-1 && j<= (gs2+gs4)-1){
          if (dy[i][j]< d/2){
             resi[i][j] = -r+(2*r*dy[i][j]*2/d);
          }else{
                 resi[i][j] = r-(2*r*(dy[i][j]-0.5*d)*2/d);
          }
          Area[i][j] = x3*x2;
          Ix[i][j] = (x2*(x3**3))/12;
          Iy[i][j] = (x3*(x2**3))/12;
        }
        area = area+ Area[i][j];
      }
    }
  }
  var ix = 0;
  var iy = 0;
  for (var i = 0; i <= x; i++){
    for (var j = 0; j <= y; j++){
      if (matrix[i][j]==1){
        iy=iy+Area[i][j]*((Math.abs(dx[i][j]-Gx))**2);
        ix=ix+Area[i][j]*((Math.abs(dy[i][j]-Gy))**2);
        resi[i][j] = resi[i][j] * 0.00125;
        
      }
    }
  }
  return [matrix, x,y,Gx, Gy,dx,dy,Area,area,Ix,Iy,ix,iy,resi];
}
function temperature(matrix,b,d,tf,gs1,gs2,gs3,dx,dy,Gx,Gy,x,y,t1,t2,t3,t4,t5,t6,t7){
var temp = [];
for(var i=0; i<=x; i++) {
    temp[i] = [];
    for(var j=0; j<=y; j++) {
      temp[i][j] = 0;
    }
}
var A1=2*(t1+t3-2*t2)/(b**2);
var B1=(t3-t1)/b;
var C1=t2;
var A2=2*(t6+t2-2*t4)/((d-2*tf)**2);
var B2=(t2-t6)/(d-2*tf);
var C2=t4;
var A3=2*(t5+t7-2*t6)/(b**2);
var B3=(t7-t5)/b;
var C3=t6;
for (var i=0; i<=x; i++){ 
    for (var j=0; j<=y; j++){
        if (matrix[i][j]==1){
            if (i<=gs1-1){
                temp[i][j]=A1*((dx[i][j]-Gx)**2)+B1*(dx[i][j]-Gx)+C1;
            }
            if (i>gs1-1 && i<=gs1+gs3-1){
                temp[i][j]=A2*((dy[i][j]-Gy)**2)+B2*(dy[i][j]-Gy)+C2;
            }
            if (i>gs1+gs3-1){
               temp[i][j]=A3*((dx[i][j]-Gx)**2)+B3*(dx[i][j]-Gx)+C3;
            }
        }
    }
}
return [temp];
}
function tempelong(temp,x,y,matrix){
var tstrain = [];
for(var i=0; i<=x; i++) {
    tstrain[i] = [];
    for(var j=0; j<=y; j++) {
      tstrain[i][j] = 0;
    }
}
for (var i=0; i<=x; i++){
    for (var j=0; j<=y; j++){
        if (matrix[i][j]==1){
            if (temp[i][j]<=750){
                tstrain[i][j]=1.2*(10**(-5))*temp[i][j]+0.4*(10**(-8))*((temp[i][j])**2)-2.416*(10**(-4));
            }
            if (temp[i][j]>750 && temp[i][j]<=860){
                tstrain[i][j]=1.1*(10**(-2));
            }
            if (temp[i][j]>860){
                tstrain[i][j]=2*(10**(-5))*temp[i][j]-6.2*(10**(-3));
            }
        }
    }
}
return [tstrain];
}
function stress(strain,i,j,temp,E,fy){
var fyt, fpt, Et,ept,eyt,c,b,a;
if (temp[i][j]<=100){
fyt = fy;
fpt = fy;
Et = E;
}
if (temp[i][j]>100 && temp[i][j]<=200){
    fyt=fy;
    fpt=fy-((temp[i][j]-100)*(fy-0.807*fy)/100);
    Et=E-((temp[i][j]-100)*(E-0.9*E)/100);
}
if (temp[i][j]>200 && temp[i][j]<=300){
    fyt=fy;
    fpt=0.807*fy-((temp[i][j]-200)*(0.807*fy-0.613*fy)/100);
    Et=0.9*E-((temp[i][j]-200)*(0.9*E-0.8*E)/100);
}
if (temp[i][j]>300 && temp[i][j]<=400){
    fyt=fy;
    fpt=0.613*fy-((temp[i][j]-300)*(0.613*fy-0.42*fy)/100);
    Et=0.8*E-((temp[i][j]-300)*(0.8*E-0.7*E)/100);
}
if (temp[i][j]>400 && temp[i][j]<=500){
    fyt=fy-((temp[i][j]-400)*(fy-0.78*fy)/100);
    fpt=0.42*fy-((temp[i][j]-400)*(0.42*fy-0.36*fy)/100);
    Et=0.7*E-((temp[i][j]-400)*(0.7*E-0.6*E)/100);
}
if (temp[i][j]>500 && temp[i][j]<=600){
    fyt=0.78*fy-((temp[i][j]-500)*(0.78*fy-0.47*fy)/100);
    fpt=0.36*fy-((temp[i][j]-500)*(0.36*fy-0.18*fy)/100);
    Et=0.6*E-((temp[i][j]-500)*(0.6*E-0.31*E)/100);
}
if (temp[i][j]>600 && temp[i][j]<=700){
    fyt=0.47*fy-((temp[i][j]-600)*(0.47*fy-0.23*fy)/100);
    fpt=0.18*fy-((temp[i][j]-600)*(0.18*fy-0.075*fy)/100);
    Et=0.31*E-((temp[i][j]-600)*(0.31*E-0.13*E)/100);
}
if (temp[i][j]>700 && temp[i][j]<=800){
    fyt=0.23*fy-((temp[i][j]-700)*(0.23*fy-0.11*fy)/100);
    fpt=0.075*fy-((temp[i][j]-700)*(0.075*fy-0.05*fy)/100);
    Et=0.13*E-((temp[i][j]-700)*(0.13*E-0.09*E)/100);
}
if (temp[i][j]>800 && temp[i][j]<=900){
    fyt=0.11*fy-((temp[i][j]-800)*(0.11*fy-0.06*fy)/100);
    fpt=0.05*fy-((temp[i][j]-800)*(0.05*fy-0.0375*fy)/100);
    Et=0.09*E-((temp[i][j]-800)*(0.09*E-0.0675*E)/100);
}
ept=fpt/Et;
eyt=0.02;
c=((fyt-fpt)**2)/((eyt-ept)*Et-2*(fyt-fpt));
b=(Math.abs(c*(eyt-ept)*Et+c**2))**0.5;
a=(Math.abs((eyt-ept)*(eyt-ept+(c/Et))))**0.5;
if (strain>=0){
    if (strain<=ept){
        f=strain*Et;
        slope=Et;
	}
    
    if (strain>ept && strain<= eyt){
        f=fpt-c+(b/a)*((Math.abs((a**2)-(eyt-strain)**2))**0.5);
        slope=(b*(eyt-strain))/(a*(Math.abs((a**2)-(eyt-strain)**2))**0.5);
	}
    if (strain>eyt){
        f=fyt;
        slope=0;
	}
}
if (strain<0){
	var negstr = -strain;
  strain = negstr;
    if (strain<=ept){
        f=strain*Et;
        slope=Et;
	}
    
    if (strain>ept && strain<= eyt){
        f=fpt-c+(b/a)*((Math.abs((a**2)-(eyt-strain)**2))**0.5);
        slope=(b*(eyt-strain))/(a*(Math.abs((a**2)-(eyt-strain)**2))**0.5);
	}
    if (strain>eyt){
        f=fyt;
        slope=0;
	}
	negstr = -f;
	f = negstr;
}
 return [f,slope];
}