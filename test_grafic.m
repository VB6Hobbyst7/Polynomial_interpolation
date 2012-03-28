function test_grafic(tip)
	switch tip
		case {0}
			# trasare grafic pentru toate cele 6 polinoame
			subplot(1,2,1);
			# trasare grafic pentru polinom lagrange, caz continuu
			eval_interpolator_c1(1,0.158,'k',0,0);
			hold on;
			# trasare grafic pentru polinom newton, caz continuu
			eval_interpolator_c1(2,0.158,'g',0,0);
			hold on;
			# trasare grafic pentru splien-uri liniare, caz continuu			
			eval_interpolator_c1(3,0.158,'b',1,0);
			hold on;
			# trasare grafic pentru spline-uri naturale, caz continuu			
			eval_interpolator_c1(4,0.158,'m',0,0);
			hold on
			# trasare grafic pentru spline-uri cubice, caz continuu			
			eval_interpolator_c1(5,0.158,'c',0,0);	
			hold on;
			# trasare grafic pentru polinom fourrier, caz continuu			
			eval_interpolator_c1(6,0.158,'y',0,0);								
			hold on;
			subplot(1,2,2);
			# trasare grafic pentru polinom lagrange, caz discret
			eval_interpolator_d1(1,10,'k',0,0);
			hold on;
			# trasare grafic pentru polinom newton, caz discret			
			eval_interpolator_d1(2,10,'g',0,0);
			hold on;
			# trasare grafic pentru splien-uri liniare, caz discret			
			eval_interpolator_d1(3,10,'b',1,0);
			hold on;
			# trasare grafic pentru spline-uri naturale,caz discret						
			eval_interpolator_d1(4,10,'m',0,0);									
			hold on;
			# trasare grafic pentru spline-uri cubice,caz discret									
			eval_interpolator_d1(5,10,'c',0,0);
			hold on;
			# trasare grafic pentru polinom fourrier,caz discret									
			eval_interpolator_d1(6,10,'y',0,0);			
			hold off;
		case {1}
			subplot(1,2,1);
			# trasare grafic pentru polinom lagrange, caz continuu			
			eval_interpolator_c1(1,0.158,'b',1,1);
			hold on;
			subplot(1,2,2);
			# trasare grafic pentru polinom lagrange, caz discret			
			eval_interpolator_d1(1,10,'b',1,1);
			hold off;
		case {2}
			subplot(1,2,1);
			# trasare grafic pentru polinom newton, caz continuu			
			eval_interpolator_c1(2,0.158,'b',1,1);
			hold on;
			subplot(1,2,2);
			# trasare grafic pentru polinom newton, caz discret						
			eval_interpolator_d1(2,10,'b',1,1);
			hold off;
		case {3}
			subplot(1,2,1);
			# trasare grafic pentru splien-uri liniare, caz continuu						
			eval_interpolator_c1(3,0.158,'b',1,1);
			hold on;
			subplot(1,2,2);
			# trasare grafic pentru splien-uri liniare, caz discret						
			eval_interpolator_d1(3,10,'b',1,1);
			hold off;
		case {4}
			subplot(1,2,1);
			# trasare grafic pentru spline-uri naturale, caz continuu						
			eval_interpolator_c1(4,0.158,'b',1,1);
			hold on;
			subplot(1,2,2);
			# trasare grafic pentru spline-uri naturale,caz discret									
			eval_interpolator_d1(4,10,'b',1,1);
			hold off;
		case {5}
			subplot(1,2,1);
			# trasare grafic pentru spline-uri cubice, caz continuu						
			eval_interpolator_c1(5,0.158,'b',1,1);
			hold on;
			subplot(1,2,2);
			# trasare grafic pentru spline-uri cubice,caz discret												
			eval_interpolator_d1(5,10,'b',1,1);
						hold off;
		case {6}
			subplot(1,2,1);
			# trasare grafic pentru polinom fourrier, caz continuu						
			eval_interpolator_c1(6,0.158,'b',1,1);
			hold on;
			subplot(1,2,2);
			# trasare grafic pentru polinom fourrier,caz discret												
			eval_interpolator_d1(6,10,'b',1,1);
						
	endswitch	
	
	
endfunction


#===============================================================================
# 		evaluare discreta
#===============================================================================

function Nk=eval_interpolator_d1(tip, ep, colr, ok, titlu)   
		maxiter=8;
		# punctele in care vor fi evaluate polinoamele de interpolare
		matr=load ("-ascii","sunspot.dat");
		x=matr(:,1);
		x=x';

		# valorile functiei f in abscisele x
		y=matr(:,2);	
		y=y';
	
		h=2*pi/300;		
		k=2;
		Nk=2^k;
		x1=y1=zeros(Nk);
		# punctele de interpolare
		index=linspace(1, 300, Nk);
		x1=x(int32(index(1:length(index))));		
		# valorile functiei f in punctele de interpolare x1
		y1=y(int32(index(1:length(index))));		
		switch   tip
			case {1}
				# caz polinom lagrange
				for l=1:300
							mypoly(l)=lagrange(Nk,x(l),x1,y1);
				endfor

			case {2}
				# caz polinom newton
				for l=1:300
					mypoly(l)=Newton(Nk,x(l),x1,y1);
				endfor
			case {3}
				# caz spline-uri liniare						
				mypoly=linear_spline(Nk,x,x1,y1);
			case {4}
				# caz spline-uri naturale							
				for l=1:300	
					mypoly(l)=ncspline( x(l),x1,y1);
				endfor
			case {5}
				# caz spline-uri cubice
				for l=1:300	
					mypoly(l)=ccspline( x(l),x1,y1);
				endfor
			case {6}
				#caz fourrier
				x1=y1=zeros(Nk);	
				index=linspace(1, 300, Nk);
				x1=x(int32(index(1:length(index))));	
				y1=y(int32(index(1:length(index))));		
				mypoly=zeros(300);
				mypoly=fourier(Nk,x,x1,y1,190);										
											
						
		endswitch	
	
		# este calculata eroarea pentru primul set de noduri de interpolare
		E1=sqrt(h*sum((y(1:300)- mypoly(1:300)).^2))/10000;	 
		#E1=sqrt(h*sum((y(1:10)- poly(1:10)).^2));	 
		while k<maxiter
			k+=1;
			Nk=2^k;
			# punctele de interpolare			
			x1=y1=zeros(Nk);	
			index=zeros(Nk);
			index=linspace(1, 300, Nk);
			x1=x(int32(index(1:length(index))));			
			# valorile functiei f in punctele de interpolare x1
			y1=y(int32(index(1:length(index))));
			
			switch   tip
					case {1}
						# caz polinom lagrange
						for l=1:300
							mypoly(l)=lagrange(Nk,x(l),x1,y1);
						endfor
					case {2}
						# caz polinom newton
						for l=1:300
						mypoly(l)=Newton(Nk,x(l),x1,y1);
						endfor
					case {3}
						# caz spline-uri liniare										
						
						mypoly=linear_spline(Nk,x,x1,y1);
					case {4}
						# caz spline-uri naturale							
							for l=1:300	
								mypoly(l)=ncspline( x(l),x1,y1);
							endfor
					case {5}
						# caz spline-uri cubice
						for l=1:300	
							mypoly(l)=ccspline( x(l),x1,y1);
						endfor
					case {6}
						#caz fourrier					
						x1=y1=zeros(Nk);	
						index=linspace(1, 300, Nk);
						x1=x(int32(index(1:length(index))));	
						y1=y(int32(index(1:length(index))));		
						mypoly=zeros(300);
						mypoly=fourier(Nk,x,x1,y1,190);						
												
						
			endswitch
   			# este calculata eroarea in nodurile de interpolare
			E2=sqrt(h*sum((y(1:300)- mypoly(1:300)).^2))/10000;	
			#E2=sqrt(h*sum((y(1:10)- poly(1:10)).^2));	
		
			
			
			# interpolantul converge daca eroare este descrescatoare 
			# si daca diferenta intre 2 erori succesive este mai mica decat un epsilon dat
			if (abs(E2 - E1) < ep &&(E2-E1)<0)
					Nk=2^k;
					break;
			endif		
			
			E1=E2;
		endwhile	
		# daca numarul de pasi atinge numarul maxim de iteratii 
		# atunci se returneaza inf(polinomul nu converge)

		
		# trasare grafic caz fourrier					
		if tip==6						
			
				if ok == 1
					plot(x1(1:Nk),y1(1:Nk),'*g;Puncte de interpolare;',x(1:300),y(1:300),'r;Functia;',1701:2000, mypoly(1:300),strcat(colr,';Fourrier;')); 
				endif
				if ok==0
					plot(1701:2000, mypoly(1:300),strcat(colr,';Fourrier;')); 
				endif
			

		else
					# legenda grafic																
				switch tip
				case {1}
					mystr=';Lagrange;';
					mystr=strcat(colr, mystr);
				case {2}
					mystr=';Newton;';
					mystr=strcat(colr, mystr);
				case {3}
					mystr=';Linear Spline;';
					mystr=strcat(colr, mystr);
				case {4}
					mystr=';Natural Spline;';
					mystr=strcat(colr, mystr);					
				case {5}
					mystr=';Cubic Spline;';
					mystr=strcat(colr, mystr);
															
				endswitch
				# trasare grafic cazuri 1-5
				if ok == 1
					plot(x1(1:Nk),y1(1:Nk),'*g;Puncte de interpolare;',1701:2000,y(1:300),'r;Functia;',1701:2000, mypoly(1:300),mystr); 
				endif
				if ok == 0
					plot(1701:2000, mypoly(1:300),mystr); 
				endif				
												
		endif
		# trasare axe de coordonate si scriere titlu
		axis([1700,2000,-200,700]);
		switch tip
				case {1}
					if Nk==256
						Nk=Inf;
					endif
					title(strcat("Evaluare discreta\nPolinom Lagrange\neps=",num2str(ep)," NK=", num2str(Nk)));
					
					
				case {2}
					if Nk==256
						Nk=Inf;
					endif
					title(strcat("Evaluare discreta\nPolinom Newton\neps=",num2str(ep)," NK=", num2str(Nk)));
				case {3}
					title(strcat("Evaluare discreta\nLinear Spline\neps=",num2str(ep)," NK=", num2str(Nk)));

				case {4}
					title(strcat("Evaluare discreta\nNatural Spline\neps=",num2str(ep)," NK=", num2str(Nk)));
				case {5}
					title(strcat("Evaluare discreta\nCubic Spline\neps=",num2str(ep)," NK=", num2str(Nk)));
				case {6}
					title(strcat("Evaluare discreta\nPolinom Fourrier\neps=",num2str(ep)," NK=", num2str(Nk)));
		endswitch
		if(titlu ==0 )title(strcat("Evaluare discreta\neps=",num2str(ep)," NK=", num2str(Nk)));
		endif
endfunction





#===============================================================================
# 		evaluare continua
#===============================================================================


function Nk=eval_interpolator_c1(tip, ep,colr, ok, titlu)   
		maxiter=11;
		# punctele in care vor fi evaluate polinoamele de interpolare
		x=linspace(-pi, pi, 1001);
		# valorile functiei f in abscisele x
		y=fct(x(1:1001));
		var=linspace(-1,1,1001);
		h=2*pi/1001;		
		k=2;
		Nk=2^k;
		x1=y1=zeros(Nk);
		# punctele de interpolare
		x1=linspace(-pi, pi,Nk);
		# valorile functiei f in punctele de interpolare x1
		y1(1:Nk)=fct(x1(1:Nk));
		switch   tip
			case {1}
				# caz polinom lagrange
				for l=1:1001
					poly(l)=lagrange(Nk,x(l),x1,y1);
				endfor
			case {2}
				# caz polinom newton
				for l=1:1001
					poly(l)=Newton(Nk,x(l),x1,y1);
				endfor
			case {3}
				# caz spline-uri liniare										
				poly=linear_spline(Nk, x,x1,y1);


			case {4}
				# caz spline-uri naturale									
				for l=1:1001	
					poly(l)=ncspline( x(l),x1,y1);
				endfor
				
			
			case {5}
				# caz spline-uri cubice													
				for l=1:1001	
					poly(l)=ccspline( x(l),x1,y1);
				endfor

			case {6}
				# caz fourrier
				x1=y1=zeros(Nk);	
				x1 = linspace(-1,1,Nk);
				y1 (1:length(x1))= fct(pi*x1(1:length(x1)));
				poly=zeros(1001);
				poly=fourier(Nk,var,x1,y1,2);
						
		endswitch	

		# este calculata eroarea pentru primul set de noduri de interpolare
		E1=sqrt(h*sum((y(1:1001)- poly(1:1001)).^2));	 
		#E1=sqrt(h*sum((y(1:10)- poly(1:10)).^2));	 
		while k<maxiter
			k+=1;
			Nk=2^k;
			# punctele de interpolare			
			x1=y1=zeros(Nk);	
			# valorile functiei f in punctele de interpolare x1					
			x1=linspace(-pi, pi,Nk );
			y1(1:Nk)=fct(x1(1:Nk));
			switch   tip
					case {1}
						# caz polinom lagrange
						for l=1:1001
							poly(l)=lagrange(Nk,x(l),x1,y1);
						endfor
					case {2}
						# caz polinom newton
						for l=1:1001
						poly(l)=Newton(Nk,x(l),x1,y1);
						endfor
					case {3}
						# caz spline-uri liniare																						
						poly=linear_spline(Nk,x,x1,y1);
						
					case {4}
						# caz spline-uri naturale														
						for l=1:1001	
							poly(l)=ncspline( x(l),x1,y1);
						endfor
					case {5}
						# caz spline-uri cubice																		
						for l=1:1001	
							poly(l)=ccspline( x(l),x1,y1);
						endfor
					case {6}
						# caz fourrier	
						x1=y1=zeros(Nk);	
						x1 = linspace(-1,1,Nk);
						y1 (1:length(x1))= fct(pi*x1(1:length(x1)));
						poly=zeros(1001);
						poly=fourier(Nk,var,x1,y1,2);
						
						
			endswitch
   			# este calculata eroarea in nodurile de interpolare
			E2=sqrt(h*sum((y(1:1001)- poly(1:1001)).^2));	
			#E2=sqrt(h*sum((y(1:10)- poly(1:10)).^2));	
			
			
			# interpolantul converge daca eroare este descrescatoare 
			# si daca diferenta intre 2 erori succesive este mai mica decat un epsilon dat
			if (abs(E2 - E1) < ep &&(E2-E1)<0)
					Nk=2^k;
					break;
			endif		
			
			E1=E2;
		endwhile
		# trasare grafic caz fourrier									
		if tip==6

				mystr=strcat(colr,';Fourrier;');
				if ok == 1
					plot(x1(1:Nk)*pi,y1(1:Nk),'*g;Punctele de interpolare;',x,y,'r;Functia;',x, poly,mystr); 
				endif	
				if ok ==0 
					plot(x, poly,mystr); 
				endif
			else
					# legenda grafic
					switch tip
					case {1}
						mystr=';Lagrange;';
						mystr=strcat(colr, mystr);
					case {2}
						mystr=';Newton;';
						mystr=strcat(colr, mystr);
					case {3}
						mystr=';Linear Spline;';
						mystr=strcat(colr, mystr);
					case {4}
						mystr=';Natural Spline;';
						mystr=strcat(colr, mystr);					
					case {5}
						mystr=';Cubic Spline;';
						mystr=strcat(colr, mystr);
															
					endswitch
					# trasare grafic cazuri 1-5										
					if ok == 1
						plot(x1(1:Nk),y1(1:Nk),'*g;Punctele de interpolare;',x,y,'r;Functia;',x, poly,mystr); 
					endif
					if ok == 0
						plot(x, poly,mystr); 
					endif												
			endif
			# trasare axe de coordonate si scriere titlu
			axis([-4,4,-0.1,0.8]);
			switch tip
				case {1}
					
					title(strcat("Evaluare continua\nPolinom Lagrange\neps=",num2str(ep)," NK=", num2str(Nk)));
				case {2}
					
					title(strcat("Evaluare continua\nPolinom Newton\neps=",num2str(ep)," NK=", num2str(Nk)));
				case {3}
					title(strcat("Evaluare continua\nLinear Spline\neps=",num2str(ep)," NK=", num2str(Nk)));
				case {4}
					title(strcat("Evaluare continua\nNatural Spline\neps=",num2str(ep)," NK=", num2str(Nk)));
				case {5}
					title(strcat("Evaluare continua\nCubic Spline\neps=",num2str(ep)," NK=", num2str(Nk)));
				case {6}
					title(strcat("Evaluare continua\nPolinom Fourrier\neps=",num2str(ep)," NK=", num2str(Nk)));
			endswitch
			if(titlu ==0 )title(strcat("Evaluare continua\neps=",num2str(ep)," NK=", num2str(Nk)));
			endif
endfunction



#===============================================================================
# 		functia pentru evaluare continua
#===============================================================================

function f=fct(x)
	f=(exp(3 * cos(x)))/(2* pi * besseli(0,3));	
endfunction


#===============================================================================
# 		Lagrange
#===============================================================================

function b = lagrange(n,a, x, y)
	
	#valoare polinom Lagrange in a
	#Intrari:
	#		a = abscisa in care se cere polinomul
	#		x = abscisele celor n+1 puncte
	#		y = ordonatele celor n+1 puncte
	#Iesiri:valoare polinom interpolare in a
	
 	b = 0;
 	for i = 1 : n
   		produs = y(i);
   		for j = 1 : n
     			if i != j
        			produs = produs * ( ( a - x( j ) ) / ( x( i ) - x( j ) ));
     			endif
   		endfor
   		b = b + produs;
 	endfor

endfunction


#===============================================================================
# 		Newton
#===============================================================================

function b=Newton(n, a, x, y)
	
	#valoare polinom Newton in a
	#Intrari:
	#		a = abscisa in care se cere polinomul
	#		x = abscisele celor n+1 puncte
	#		y = ordonatele celor n+1 puncte
	#Iesiri:valoare polinom interpolare in a
	
	for k=1:n-1
		y(k+1:n)=(y(k+1:n)-y(k))./(x(k+1:n)-x(k));
	endfor
	c=y(:);
	 b=c(1);
	 p=1;
	for i=2:n
		p=(a-x(i-1)).*p;
		b=b+p*c(i);
	endfor

endfunction


#===============================================================================
# 		Fourrier
#===============================================================================

function  yInt=fourier(n,xInt,x,y,L)	
	n = length(x);
	m = n/2; 
	xx = [x(m+1:n),x(1:m)]';	
	yy = [y(m+1:n),y(1:m)]';	
	for j = 0 : m
		a(j+1) = 2*yy'*cos(2*pi*j*xx/L)/n;
		b(j+1) = 2*yy'*sin(2*pi*j*xx/L)/n;
	endfor
	a';
	b'; 
	yInt = 0.5*a(1)*ones(1,length(xInt));
	for j = 1 : (m-1)
		yInt = yInt + a(1+j)*cos(2*pi*j*xInt/L) + b(1+j)*sin(2*pi*j*xInt/L);
	endfor
	yInt = yInt + a(m+1)*cos(2*pi*m*xInt/L);
endfunction



#===============================================================================
# 		Linear Spline
#===============================================================================

function f=linear_spline(n,xInt,x,y)

	n = length(x)-1; 
	a = y(1:n);
	b = (y(2:n+1)-y(1:n))./(x(2:n+1)-x(1:n));	 
	for j = 1 : length(xInt)
	  if xInt(j) ~= x(n+1)
			iInt(j) = sum(x <= xInt(j));
	  else
	  iInt(j) = n;
	  end
	end
	yInt = a(iInt) + b(iInt).*(xInt-x(iInt));
	f=yInt;
endfunction



#===============================================================================
# 		Natural cubic spline
#===============================================================================
function yy=ncspline(xx, x,y)
	h=diff(x);
	n=length(x)-1;  
	a(1:n+1)=y(1:n+1);
	A=sparse(2:n,1:n-1,h(1:n-1),n+1,n+1) + ...
	  sparse(2:n,3:n+1,h(2:n),n+1,n+1) + ...
	  sparse(2:n,2:n,2*(h(1:n-1)+h(2:n)),n+1,n+1);
	A(1,1)=1; A(n+1,n+1)=1;
	b=[0,3./h(2:n).*(a(3:n+1)-a(2:n))-3./h(1:n-1).*(a(2:n)-a(1:n-1)),0]';
	c=(A\b)';
	b=(a(2:n+1)-a(1:n))./h-h./3.*(2*c(1:n)+c(2:n+1));
	d=(c(2:n+1)-c(1:n))./(3*h);
	c=c(1:n);
	for i=1:n
		if xx>=x(i) & xx<=x(i+1);
			  yy=a(i)+b(i)*(xx-x(i))+c(i)*(xx-x(i)).^2+d(i)*(xx-x(i)).^3;
			  break;
		endif
	endfor
endfunction



#===============================================================================
# 		Cubic spline
#===============================================================================
function yy=ccspline(xx,x,y)
		h=diff(x);
		n=length(x)-1;
		fa=(y(2)-y(1))/(x(2)-x(1));
		fb=(y(n+1)-y(n))/(x(n+1)-x(n));
		a(1:n+1)=y(1:n+1);
		A=sparse(2:n+1,1:n,h,n+1,n+1) + ...
		  sparse(1:n,2:n+1,h,n+1,n+1) + ...
		  sparse(2:n,2:n,2*(h(1:n-1)+h(2:n)),n+1,n+1);
		A(1,1)=2*h(1); A(n+1,n+1)=2*h(n);
		b=[3./h(1)*(a(2)-a(1))-3*fa,3./h(2:n).*(a(3:n+1)-a(2:n))-3./h(1:n-1).*(a(2:n)-a(1:n-1)),3*fb-3/h(n)*(a(n+1)-a(n))]';
		c=(A\b)';
		b=(a(2:n+1)-a(1:n))./h-h./3.*(2*c(1:n)+c(2:n+1));
		d=(c(2:n+1)-c(1:n))./(3*h);
		c=c(1:n);
		for i=1:n
		if xx>=x(i) & xx<=x(i+1);
			  yy=a(i)+b(i)*(xx-x(i))+c(i)*(xx-x(i)).^2+d(i)*(xx-x(i)).^3;
			  break;
		endif
	endfor
endfunction

