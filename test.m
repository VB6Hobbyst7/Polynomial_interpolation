function test()
	# matrice corespunzatoare numarului de noduri de interpolare 
	# pentru care fiecare interpolant converge
	m=cell(2, 6);
	m(1,1)=eval_interpolator_c(1,0.158);
	m(1,2)=eval_interpolator_c(2,0.158);
	m(1,3)=eval_interpolator_c(3,0.158);
	m(1,4)=eval_interpolator_c(4,0.158);
	m(1,5)=eval_interpolator_c(5,0.158);
	m(1,6)=eval_interpolator_c(6,0.158);
	m(2,1)=eval_interpolator_d(1,10);
	m(2,2)=eval_interpolator_d(2,10);
	m(2,3)=eval_interpolator_d(3,10);
	m(2,4)=eval_interpolator_d(4,10);
	m(2,5)=eval_interpolator_d(5,10);
	m(2,6)=eval_interpolator_d(6,10);
	# afisare matrice
	disp(m);

endfunction
