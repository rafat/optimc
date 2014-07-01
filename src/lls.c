/*
 * lls.c
 *
 *  Created on: Apr 14, 2014
 *      Author: Rafat Hussain
 */

#include "lls.h"

int lls_normal(double *A,double *b,int M,int N,double *x) {
	int retcode,sc,i,j,c1,l;
	double sum;
	double *AT,*d,*C,*y,*CT;
	// M - data points
	// N - Number of parameters
	// A - MXN; b size - M vector; AT - NXM

	AT = (double*) malloc(sizeof(double) * M * N);
	d = (double*) malloc(sizeof(double) * N);
	C = (double*) malloc(sizeof(double) * N * N);
	y = (double*) malloc(sizeof(double) * N);
	CT = (double*) malloc(sizeof(double) * N * N);

	retcode = 0;

	mtranspose(A,M,N,AT);
	mmult(AT,b,d,N,M,1);
	mmult(AT,A,C,N,M,N);

	sc = chol(C,N);
	if (sc == -1) {
		return -1;
	}

	mtranspose(C,N,N,CT);
	//Forward Substitution

	y[0] = d[0]/CT[0];
	for(i = 1; i < N; ++i) {
		sum = 0.;
		c1 = i*N;
		for(j = 0; j < i; ++j) {
			sum += y[j] * CT[c1 + j];
		}
		y[i] = (d[i] - sum)/CT[c1+i];
	}

	//Back Substitution

	x[N - 1] = y[N - 1]/C[N * N - 1];

	for (i = N - 2; i >= 0; i--) {
		sum = 0.;
		c1 = i*(N+1);
		l=0;
		for(j = i+1; j < N;j++) {
			l++;
			sum += C[c1 + l] * x[j];
		}
		x[i] = (y[i] - sum) / C[c1];
	}

	free(AT);
	free(d);
	free(C);
	free(y);
	return retcode;
}

int lls_qr(double *Ai,double *bi,int M,int N,double *xo) {
	int j,i,k,u,t,retcode,c1,l;
	double *x,*v,*AT,*w,*bvec,*b,*A,*R;
	double beta,sum;

	retcode = 0;

	if (M < N) {
			printf("M should be greater than or equal to N");
			exit(-1);
	}
	x = (double*) malloc(sizeof(double) * M);
	b = (double*) malloc(sizeof(double) * M);
	bvec = (double*) malloc(sizeof(double) * N);
	v = (double*) malloc(sizeof(double) * M);
	AT = (double*) malloc(sizeof(double) * M * N);
	A = (double*) malloc(sizeof(double) * M * N);
	w = (double*) malloc(sizeof(double) * 1);
	R = (double*) malloc(sizeof(double) * N * N);

	for(j = 0; j < M;++j) {
		b[j] = bi[j];
	}
	for(j = 0; j < M*N;++j) {
		A[j] = Ai[j];
	}

	for(j = 0; j < N;++j) {
		for(i=j;i < M;++i) {
			x[i-j] = A[i*N+j];

		}

		beta = house(x,M-j,v);
		bvec[j] = beta;

		for (i=j; i < M; i++) {
			t = i * N;
			u = 0;
			for (k=j; k < N; k++) {
				AT[u+i-j] = A[k+t];
				u+=(M-j);

			}

		}


		mmult(AT,v,w,N-j,M-j,1);
		scale(w,N-j,1,beta);
		mmult(v,w,AT,M-j,1,N-j);
		for (i=j; i < M; i++) {
			t = i *N;
			for (k=j; k < N; k++) {
				A[t+k] -= AT[(i-j)*(N-j) + k - j];
			}
		}
		if (j < M) {

			for(i=j+1;i < M;++i) {
				A[i*N+j] = v[i-j];
			}
		}

	}

	for(i = 0; i < N;++i) {
		t = i *N;
		for(j = 0; j < N;++j) {
			if (i > j) {
				R[t+j] = 0.;
			} else {
				R[t+j] = A[t+j];
			}
		}
	}


	for(j = 0; j < N;++j) {
		v[j] = 1;
		for(i = j+1; i < M;++i) {
			v[i] = A[i * N + j];//edit
		}
		mmult(b+j,v+j,w,1,M-j,1);
		*w = *w * bvec[j];
		for(i = j; i < M;++i) {
			v[i] = *w * v[i];
		}
		for(i = j; i < M;++i) {
			b[i] = b[i] - v[i];
		}


	}

	//mdisplay(b,1,M);

	//back substitution

	xo[N - 1] = b[N - 1]/R[N * N - 1];

	for (i = N - 2; i >= 0; i--) {
		sum = 0.;
		c1 = i*(N+1);
		l=0;
		for(j = i+1; j < N;j++) {
			l++;
			sum += R[c1 + l] * xo[j];
		}
		xo[i] = (b[i] - sum) / R[c1];
	}

	free(x);
	free(v);
	free(AT);
	free(w);
	free(bvec);
	free(R);
	free(b);
	free(A);

	return retcode;
}

void bidiag(double *A, int M, int N) {
	int j,i,k,u,t;
	double *x,*v,*AT,*w;
	double beta;

	if (M < N) {
			printf("M should be greater than or equal to N");
			exit(1);
	}
	x = (double*) malloc(sizeof(double) * M);
	v = (double*) malloc(sizeof(double) * M);
	AT = (double*) malloc(sizeof(double) * M * N);
	w = (double*) malloc(sizeof(double) * M * M);


	for(j = 0; j < N;++j) {
		for(i=j;i < M;++i) {
			x[i-j] = A[i*N+j];

		}

		beta = house(x,M-j,v);

		for (i=j; i < M; i++) {
			t = i * N;
			u = 0;
			for (k=j; k < N; k++) {
				AT[u+i-j] = A[k+t];
				u+=(M-j);

			}

		}


		mmult(AT,v,w,N-j,M-j,1);
		scale(w,N-j,1,beta);
		mmult(v,w,AT,M-j,1,N-j);
		for (i=j; i < M; i++) {
			t = i *N;
			for (k=j; k < N; k++) {
				A[t+k] -= AT[(i-j)*(N-j) + k - j];
			}
		}

		for(i=j+1;i < M;++i) {
			A[i*N+j] = v[i-j];
		}

		if (j < N - 2) {
			for(i=j+1;i < N;++i) {
				x[i-j-1] = A[j*N+i];
			}
			beta = house(x,N-j-1,v);

			for (i=j; i < M; i++) {
				t = i * N;
				u = (i-j) *(N-j-1);
				for (k=j+1; k < N; k++) {
					AT[u+k-j-1] = A[k+t];
				}
			}

			mmult(AT,v,w,M-j,N-j-1,1);
			scale(w,M-j,1,beta);
			mmult(w,v,AT,M-j,1,N-j-1);

			for (i=j; i < M; i++) {
				t = i * N;
				u = (i-j) *(N-j-1);
				for (k=j+1; k < N; k++) {
					A[k+t] -= AT[u+k-j-1];
				}
			}
			u = 1;
			t = j*N;
			for(i = j + 2; i < N;++i) {
				A[t+i] = v[u];
				u++;
			}

		}

	}


	free(x);
	free(v);
	free(AT);

}

void bidiag_orth(double *A, int M, int N,double *U,double *V) {
	int j,i,k,u,t;
	double *x,*v,*AT,*w;
	double beta;

	if (M < N) {
			printf("M should be greater than or equal to N");
			exit(1);
	}
	x = (double*) malloc(sizeof(double) * M);
	v = (double*) malloc(sizeof(double) * M);
	AT = (double*) malloc(sizeof(double) * M * M);
	w = (double*) malloc(sizeof(double) * M * M);


	eye(U,M);
	eye(V,N);

	for(j = 0; j < N;++j) {
		for(i=j;i < M;++i) {
			x[i-j] = A[i*N+j];

		}

		beta = house(x,M-j,v);
		//mdisplay(v,M-j,1);

		for (i=j; i < M; i++) {
			t = i * N;
			u = 0;
			for (k=j; k < N; k++) {
				AT[u+i-j] = A[k+t];
				u+=(M-j);

			}

		}


		mmult(AT,v,w,N-j,M-j,1);
		scale(w,N-j,1,beta);
		mmult(v,w,AT,M-j,1,N-j);
		for (i=j; i < M; i++) {
			t = i *N;
			for (k=j; k < N; k++) {
				A[t+k] -= AT[(i-j)*(N-j) + k - j];
			}
		}

		for(i=j+1;i < M;++i) {
			A[i*N+j] = v[i-j];
		}


		for (i=j; i < M; i++) {
			t = i * M;
			u = (i-j) *(M-j);
			for (k=j; k < M; k++) {
				AT[u+k-j] = U[k+t];
			}
		}

		mmult(AT,v,w,M-j,M-j,1);
		scale(w,M-j,1,beta);
		mmult(w,v,AT,M-j,1,M-j);

		for (i=j; i < M; i++) {
			t = i * M;
			u = (i-j) *(M-j);
			for (k=j; k < M; k++) {
				U[k+t] -= AT[u+k-j];
			}
		}

		//mdisplay(U,M,M);

		if (j < N - 2) {
			for(i=j+1;i < N;++i) {
				x[i-j-1] = A[j*N+i];
			}
			beta = house(x,N-j-1,v);

			for (i=j; i < M; i++) {
				t = i * N;
				u = (i-j) *(N-j-1);
				for (k=j+1; k < N; k++) {
					AT[u+k-j-1] = A[k+t];
				}
			}

			mmult(AT,v,w,M-j,N-j-1,1);
			scale(w,M-j,1,beta);
			mmult(w,v,AT,M-j,1,N-j-1);

			for (i=j; i < M; i++) {
				t = i * N;
				u = (i-j) *(N-j-1);
				for (k=j+1; k < N; k++) {
					A[k+t] -= AT[u+k-j-1];
				}
			}
			u = 1;
			t = j*N;
			for(i = j + 2; i < N;++i) {
				A[t+i] = v[u];
				u++;
			}

		}

	}


	free(x);
	free(v);
	free(AT);

}


int svd_gr(double *A,int M,int N,double *U,double *V,double *q) {
	int i,j,k,l,t,t2,ierr,cancel,iter,l1;
	double eps,g,x,s,temp,f,h,scale,c,y,z;
	double *e;
	/*
     THIS SUBROUTINE IS THE MODIFIED C TRANSLATION OF THE
     EISPACK FORTRAN TRANSLATION OF THE ALGOL PROCEDURE SVD,
     NUM. MATH. 14, 403-420(1970) BY GOLUB AND REINSCH.
     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
	 */
	/*
	 * U = MXN
	 * V - NXN
	 * Q - NX1
	 */

	e = (double*) malloc(sizeof(double) * N);
	ierr = 0;
	eps = macheps();
	g = scale = x = 0.0;

	for(i = 0; i < M*N;++i) {
		U[i] = A[i];
	}

	for(i = 0; i < N;++i) {
		l = i+1;
		e[i] = scale * g;
		g = 0.0;
		s = 0.0;
		scale = 0.0;

		if (i < M) {
			for(k = i; k < M;++k) {
				scale += fabs(U[k*N+i]);
			}

			if (scale != 0.0) {
				for(k = i; k < M;++k) {
					t = k * N;
					U[t+i] /= scale;
					temp = U[t+i];
					s += temp*temp;
				}
				f = U[i*N+i];
				g = (f < 0) ? sqrt(s) : -sqrt(s);
				h = f * g - s;
				U[i*N+i] = f - g;

				if (i < N - 1) {
					for(j = l; j < N;++j) {
						s = 0.0;
						for(k = i; k < M;++k) {
							t = k * N;
							s += U[t+i]*U[t+j];
						}
						f = s / h;
						for(k = i; k < M;++k) {
							t = k * N;
							U[t+j] += f * U[t+i];
						}
					}
				}
				for(k = i; k < M;++k) {
					t = k * N;
					U[t+i] *= scale;
				}
			}
		}
        q[i] = scale * g;
        g = 0.0;
        s = 0.0;
        scale = 0.0;

        if (i < M && i != N - 1) {
        	t = i *N;
        	for(k = l; k < M;++k) {
        		scale += fabs(U[t+k]);
        	}
        	if (scale != 0.0) {
        		for(k = l; k < N;++k) {
        			U[t+k] /= scale;
        			temp = U[t+k];
        			s = s + temp*temp;
        		}
        		f = U[t+l];
        		g = (f < 0) ? sqrt(s) : -sqrt(s);
                h = f * g - s;
                U[t+l] = f - g;
                for(k = l;k < N;++k) {
                	e[k] = U[t+k] / h;
                }

				for (j = l; j < M; j++) {
					s = 0.0;
					t2 = j * N;
					for (k = l; k < N; k++) {
						s += U[t2+k] * U[t+k];
					}
					for (k = l; k < N; k++) {
						U[t2+k] += s * e[k];
					}
				}
                for (k = l; k < N; k++)
                    U[t+k] *= scale;
        	}

        }

        temp = fabs(q[i]) + fabs(e[i]);

        if (x < temp) {
        	x = temp;
        }
	}

	//Accumulating Right Hand Transformations

		for(i = N - 1;i >= 0;--i) {
			t = i * N;
			if (i < N - 1) {
				if (g != 0.0) {
					h = U[t+i+1] * g;
					for(j = l;j < N;++j) {
						V[j*N+i] = U[t+j] / h;
					}
					for(j = l;j < N;++j) {
						s = 0.0;
						for(k = l; k < N;++k) {
							s += U[t+k] * V[k*N+j];
						}
						for(k = l; k < N;++k) {
							V[k*N+j] += (s * V[k*N+i]);
						}
					}
				}
				for(j = l; j < N;++j) {
					V[t+j] = V[j*N+i] = 0.0;
				}
			}
		    V[t+i] = 1.0;
			g = e[i];
			l = i;
		}



	//Accumulating Left Hand Transformations

		for(i = N - 1;i >= 0;--i) {
			t = i * N;
			l = i+1;
			g = q[i];

			if (i < N - 1) {
				for(j = l;j < N;++j) {
					U[t+j] = 0.0;
				}
			}

			if (g != 0.0) {
				if (i != N - 1) {
					//h = U[t+i] * g;
					for(j = l;j < N;++j) {
						s = 0.0;
						for(k = l; k < M;++k) {
							s += (U[k*N+i] * U[k*N+j]);
						}
						f = (s / U[t+i]) / g;
						for(k = i; k < M;++k) {
							U[k*N+j] += (f * U[k*N+i]);
						}
					}
				}
				for(j = i; j < M;++j) {
					U[j*N+i] = U[j*N+i] / g;
				}
			} else {
				for(j = i; j < M;++j) {
					U[j*N+i] = 0.0;
				}
			}

			U[t+i] += 1.0;
		}
	//	mdisplay(U,M,N);

		eps = eps * x;

		for(k = N - 1; k >= 0; --k) {
			iter = 0;

			while(1) {
				iter++;
				if (iter > SVDMAXITER) {
					printf("Convergence Not Achieved \n");
					return 15;
				}

				cancel = 1;
				for(l = k; l >= 0; --l) {
					if (fabs(e[l]) <= eps) {
						cancel = 0; //test f convergence
						break;
					}
					if (fabs(q[l-1]) <= eps) {
						//Cancel
						break;
					}
				}
				if (cancel) {
					c = 0.0;
					s = 1.0;
					l1 = l - 1;
					for(i = l; i <= k;++i) {
						f = s*e[i];
						e[i] *= c;
						if (fabs(f) <= eps) {
							break;
						}
						g = q[i];
						h = q[i] = hypot(f,g);
						c = g/h;
						s = -f/h;
						for(j = 0; j < M;++j) {
							t = j * N;
							y = U[t+l1];
							z = U[t+i];

							U[t+l1] = y * c + z * s;
							U[t+i] = z * c - y * s;
						}
					}
				}
				z = q[k];
				if (l != k) {
					x = q[l];
					y = q[k-1];
					g = e[k-1];
					h = e[k];
					f = 0.5 * (((g + z) / h) * ((g - z) / y) + y / h - h / y);
					g = hypot(f,1.0);
					if (f < 0.0) {
						temp = f - g;
					} else {
						temp = f+g;
					}
					f = x - (z / x) * z + (h / x) * (y / temp - h);

					//Next QR Transformation

					c = s = 1.0;
					for(i = l+1; i <= k;++i) {
						g = e[i];
						y = q[i];
						h = s * g;
						g = c * g;
						e[i-1] = z = hypot(f,h);
	                    c = f / z;
	                    s = h / z;
	                    f = x * c + g * s;
	                    g = g * c - x * s;
	                    h = y * s;
	                    y *= c;
	                    for(j = 0; j < N;++j) {
	                    	t = j * N;
	                        x = V[t+i-1];
	                        z = V[t+i];
	                        V[t+i-1] = x * c + z * s;
	                        V[t+i] = z * c - x * s;
	                    }
	                    q[i-1] = z = hypot(f,h);
	                    if (z != 0.0) {
	                        c = f / z;
	                        s = h / z;
	                    }
	                    f = c * g + s * y;
	                    x = c * y - s * g;
	                    for(j = 0; j < M;++j) {
	                    	t = j * N;
	                        y = U[t+i-1];
	                        z = U[t+i];
	                        U[t+i-1] = y * c + z * s;
	                        U[t+i] = z * c - y * s;
	                    }
					}
	                    e[l] = 0.0;
	                    e[k] = f;
	                    q[k] = x;

				} else {
					//convergence
	                if (z < 0.0) {
	                    q[k] = -z;
	                    for (j = 0; j < N; j++) {
	                    	t = j *N;
	                        V[t+k] = -V[t+k];
	                    }
	                }
	                break;
				}
			}
		}

		svd_sort(U,M,N,V,q);

		free(e);
		return ierr;
}


int svd_gr2(double *A,int M,int N,double *U,double *V,double *q) {
	int i,j,k,l,t,t2,ierr,cancel,iter,l1;
	double eps,g,x,s,temp,f,h,tol,c,y,z;
	double *e;
	/*
     THIS SUBROUTINE IS THE MODIFIED C TRANSLATION OF THE
     EISPACK FORTRAN TRANSLATION OF THE ALGOL PROCEDURE SVD,
     NUM. MATH. 14, 403-420(1970) BY GOLUB AND REINSCH.
     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
	 */
	/*
	 * U = MXN
	 * V - NXN
	 * Q - NX1
	 */
	if (M < N) {
		printf("Rows (M) should be greater than Columns (B) \n");
		printf("Retry By Transposing the Input Matrix");
		return -1;
	}
	e = (double*) malloc(sizeof(double) * N);
	ierr = 0;
	eps = macheps();
	tol = eps * eps;
	g = x = 0.0;

	for(i = 0; i < M*N;++i) {
		U[i] = A[i];
	}

	for(i = 0; i < N;++i) {
		l = i+1;
		e[i] = g;
		s = 0.0;

		for(k = i; k < M;++k) {
			t = k * N;
			temp = U[t+i];
			s += temp*temp;
		}
		if (s < tol) {
			g = 0.0;
		} else {
			f = U[i*N+i];
			g = (f < 0) ? sqrt(s) : -sqrt(s);
			h = f * g - s;
			U[i*N+i] = f - g;

			for(j = l; j < N;++j) {
				s = 0.0;
				for(k = i; k < M;++k) {
					t = k * N;
					s += (U[t+i]*U[t+j]);
				}
				f = s / h;
				for(k = i; k < M;++k) {
					t = k * N;
					U[t+j] += (f * U[t+i]);
				}
			}

		}

        q[i] = g;
        s = 0.0;
        t = i * N;
    	for(k = l; k < N;++k) {
    		temp = U[t+k];
    		s = s + temp*temp;
    	}
        if (s < tol) {
        	g = 0.0;
        } else {
        	f = U[t+l];
			g = (f < 0) ? sqrt(s) : -sqrt(s);
			h = f * g - s;
            U[t+l] = f - g;
            for(k = l;k < N;++k) {
              	e[k] = U[t+k] / h;
            }

            for (j = l; j < M; j++) {
                s = 0.0;
                t2 = j * N;
                for (k = l; k < N; k++) {
                     s += U[t2+k] * U[t+k];
                }
                for (k = l; k < N; k++) {
                     U[t2+k] += s * e[k];
                }
            }

        }

        temp = fabs(q[i]) + fabs(e[i]);

        if (x < temp) {
        	x = temp;
        }
	}


//Accumulating Right Hand Transformations

	for(i = N - 1;i >= 0;--i) {
		t = i * N;
		if (i < N - 1) {
			if (g != 0.0) {
				h = U[t+i+1] * g;
				for(j = l;j < N;++j) {
					V[j*N+i] = U[t+j] / h;
				}
				for(j = l;j < N;++j) {
					s = 0.0;
					for(k = l; k < N;++k) {
						s += U[t+k] * V[k*N+j];
					}
					for(k = l; k < N;++k) {
						V[k*N+j] += (s * V[k*N+i]);
					}
				}
			}
			for(j = l; j < N;++j) {
				V[t+j] = V[j*N+i] = 0.0;
			}
		}
	    V[t+i] = 1.0;
		g = e[i];
		l = i;
	}



//Accumulating Left Hand Transformations

	for(i = N - 1;i >= 0;--i) {
		t = i * N;
		l = i+1;
		g = q[i];

		if (i < N - 1) {
			for(j = l;j < N;++j) {
				U[t+j] = 0.0;
			}
		}

		if (g != 0.0) {
			if (i != N - 1) {
				//h = U[t+i] * g;
				for(j = l;j < N;++j) {
					s = 0.0;
					for(k = l; k < M;++k) {
						s += (U[k*N+i] * U[k*N+j]);
					}
					f = (s / U[t+i]) / g;
					for(k = i; k < M;++k) {
						U[k*N+j] += (f * U[k*N+i]);
					}
				}
			}
			for(j = i; j < M;++j) {
				U[j*N+i] = U[j*N+i] / g;
			}
		} else {
			for(j = i; j < M;++j) {
				U[j*N+i] = 0.0;
			}
		}

		U[t+i] += 1.0;
	}
//	mdisplay(U,M,N);

	eps = eps * x;

	for(k = N - 1; k >= 0; --k) {
		iter = 0;

		while(1) {
			iter++;
			if (iter > SVDMAXITER) {
				printf("Convergence Not Achieved \n");
				return 15;
			}

			cancel = 1;
			for(l = k; l >= 0; --l) {
				if (fabs(e[l]) <= eps) {
					cancel = 0; //test f convergence
					break;
				}
				if (fabs(q[l-1]) <= eps) {
					//Cancel
					break;
				}
			}
			if (cancel) {
				c = 0.0;
				s = 1.0;
				l1 = l - 1;
				for(i = l; i <= k;++i) {
					f = s*e[i];
					e[i] *= c;
					if (fabs(f) <= eps) {
						break;
					}
					g = q[i];
					h = q[i] = hypot(f,g);
					c = g/h;
					s = -f/h;
					for(j = 0; j < M;++j) {
						t = j * N;
						y = U[t+l1];
						z = U[t+i];

						U[t+l1] = y * c + z * s;
						U[t+i] = z * c - y * s;
					}
				}
			}
			z = q[k];
			if (l != k) {
				x = q[l];
				y = q[k-1];
				g = e[k-1];
				h = e[k];
				f = 0.5 * (((g + z) / h) * ((g - z) / y) + y / h - h / y);
				g = hypot(f,1.0);
				if (f < 0.0) {
					temp = f - g;
				} else {
					temp = f+g;
				}
				f = x - (z / x) * z + (h / x) * (y / temp - h);

				//Next QR Transformation

				c = s = 1.0;
				for(i = l+1; i <= k;++i) {
					g = e[i];
					y = q[i];
					h = s * g;
					g = c * g;
					e[i-1] = z = hypot(f,h);
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = g * c - x * s;
                    h = y * s;
                    y *= c;
                    for(j = 0; j < N;++j) {
                    	t = j * N;
                        x = V[t+i-1];
                        z = V[t+i];
                        V[t+i-1] = x * c + z * s;
                        V[t+i] = z * c - x * s;
                    }
                    q[i-1] = z = hypot(f,h);
                    if (z != 0.0) {
                        c = f / z;
                        s = h / z;
                    }
                    f = c * g + s * y;
                    x = c * y - s * g;
                    for(j = 0; j < M;++j) {
                    	t = j * N;
                        y = U[t+i-1];
                        z = U[t+i];
                        U[t+i-1] = y * c + z * s;
                        U[t+i] = z * c - y * s;
                    }
				}
                    e[l] = 0.0;
                    e[k] = f;
                    q[k] = x;

			} else {
				//convergence
                if (z < 0.0) {
                    q[k] = -z;
                    for (j = 0; j < N; j++) {
                    	t = j *N;
                        V[t+k] = -V[t+k];
                    }
                }
                break;
			}
		}
	}

	svd_sort(U,M,N,V,q);

	free(e);
	return ierr;
}

int minfit(double *AB,int M,int N,int P,double *q) {
	int i,j,k,l,t,t2,ierr,cancel,iter,l1,np,n1;
	double eps,g,x,s,temp,f,h,tol,c,y,z;
	double *e;
	/*
     THIS SUBROUTINE IS THE MODIFIED C TRANSLATION OF THE
     EISPACK FORTRAN TRANSLATION OF THE ALGOL PROCEDURE MINFIT,
     NUM. MATH. 14, 403-420(1970) BY GOLUB AND REINSCH.
     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
	 */
	/*
	 * AB MAX(M,N) X (N+P)
	 * q N X 1
	 *
	 * On Input AB = [A:B]
	 * On Output AB = [V:C] where C = U'*B
	 */

	e = (double*) malloc(sizeof(double) * N);
	ierr = 0;
	eps = macheps();
	tol = eps * eps;
	g = x = 0.0;
	np = N + P;

	for(i = 0; i < N;++i) {
		l = i+1;
		e[i] = g;
		s = 0.0;

		for(k = i; k < M;++k) {
			t = k * np;
			temp = AB[t+i];
			s += temp*temp;
		}
		if (s < tol) {
			g = 0.0;
		} else {
			f = AB[i*np+i];
			g = (f < 0) ? sqrt(s) : -sqrt(s);
			h = f * g - s;
			AB[i*np+i] = f - g;

			for(j = l; j < np;++j) {
				s = 0.0;
				for(k = i; k < M;++k) {
					t = k * np;
					s += (AB[t+i]*AB[t+j]);
				}
				f = s / h;
				for(k = i; k < M;++k) {
					t = k * np;
					AB[t+j] += (f * AB[t+i]);
				}
			}

		}

        q[i] = g;
        s = 0.0;
        t = i * np;
        if (i < M) {
			for(k = l; k < N;++k) {
				temp = AB[t+k];
				s = s + temp*temp;
			}
        }
        if (s < tol) {
        	g = 0.0;
        } else {
        	f = AB[t+i+1];
			g = (f < 0) ? sqrt(s) : -sqrt(s);
			h = f * g - s;
            AB[t+i+1] = f - g;
            for(k = l;k < N;++k) {
              	e[k] = AB[t+k] / h;
            }

            for (j = l; j < M; j++) {
                s = 0.0;
                t2 = j * np;
                for (k = l; k < N; k++) {
                     s += AB[t2+k] * AB[t+k];
                }
                for (k = l; k < N; k++) {
                     AB[t2+k] += s * e[k];
                }
            }

        }

        temp = fabs(q[i]) + fabs(e[i]);

        if (x < temp) {
        	x = temp;
        }
	}

	//Accumulating Right Hand Transformations

	for(i = N - 1;i >= 0;--i) {
		t = i * np;
		if (g != 0.0) {
			h = AB[t+i+1] * g;
			for(j = l;j < N;++j) {
				AB[j*np+i] = AB[t+j] / h;
			}
			for(j = l;j < N;++j) {
				s = 0.0;
				for(k = l; k < N;++k) {
					s += AB[t+k] * AB[k*np+j];
				}
				for(k = l; k < N;++k) {
					AB[k*np+j] += (s * AB[k*np+i]);
				}
			}
		}
		for(j = l; j < N;++j) {
			AB[t+j] = AB[j*np+i] = 0.0;
		}

		AB[t+i] = 1.0;
		g = e[i];
		l = i;
	}

	eps = eps *x;
	n1 = N + 1;

	for(i = M; i < N;++i) {
		for(j = n1 - 1; j < np;++j)  {
			AB[i * np +j] = 0.0;
		}
	}

	for(k = N - 1; k >= 0; --k) {
		iter = 0;

		while(1) {
			iter++;
			if (iter > SVDMAXITER) {
				printf("Convergence Not Achieved \n");
				return 15;
			}

			cancel = 1;
			for(l = k; l >= 0; --l) {
				if (fabs(e[l]) <= eps) {
					cancel = 0; //test f convergence
					break;
				}
				if (fabs(q[l-1]) <= eps) {
					//Cancel
					break;
				}
			}
			if (cancel) {
				c = 0.0;
				s = 1.0;
				l1 = l - 1;
				for(i = l; i <= k;++i) {
					f = s*e[i];
					e[i] *= c;
					if (fabs(f) <= eps) {
						break;
					}
					g = q[i];
					h = q[i] = hypot(f,g);
					c = g/h;
					s = -f/h;
					for(j = n1-1; j < np;++j) {
						y = AB[l1 * np + j];
						z = AB[i * np + j];

						AB[l1 * np + j] = y * c + z * s;
						AB[i * np + j] = z * c - y * s;
					}
				}
			}
			z = q[k];
			if (l != k) {
				x = q[l];
				y = q[k-1];
				g = e[k-1];
				h = e[k];
				f = 0.5 * (((g + z) / h) * ((g - z) / y) + y / h - h / y);
				g = hypot(f,1.0);
				if (f < 0.0) {
					temp = f - g;
				} else {
					temp = f+g;
				}
				f = x - (z / x) * z + (h / x) * (y / temp - h);

				//Next QR Transformation

				c = s = 1.0;
				for(i = l+1; i <= k;++i) {
					g = e[i];
					y = q[i];
					h = s * g;
					g = c * g;
					e[i-1] = z = hypot(f,h);
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = g * c - x * s;
                    h = y * s;
                    y *= c;
                    for(j = 0; j < N;++j) {
                    	t = j * np;
                        x = AB[t+i-1];
                        z = AB[t+i];
                        AB[t+i-1] = x * c + z * s;
                        AB[t+i] = z * c - x * s;
                    }
                    q[i-1] = z = hypot(f,h);
                    if (z != 0.0) {
                        c = f / z;
                        s = h / z;
                    }
                    f = c * g + s * y;
                    x = c * y - s * g;
                    for(j = n1-1; j < np;++j) {
                        y = AB[(i-1) * np + j];
                        z = AB[i * np + j];
                        AB[(i-1) * np + j] = y * c + z * s;
                        AB[i * np + j] = z * c - y * s;
                    }
				}
                    e[l] = 0.0;
                    e[k] = f;
                    q[k] = x;

			} else {
				//convergence
                if (z < 0.0) {
                    q[k] = -z;
                    for (j = 0; j < N; j++) {
                    	t = j * np;
                        AB[t+k] = -AB[t+k];
                    }
                }
                break;
			}
		}
	}

	free(e);
	return ierr;

}

int lls_svd2(double *Ai,double *bi,int M,int N,double *xo) {
	int retcode,P,mnmax,np,i,j,t,t2;
	double *AB,*q,*V,*C;
	double eps;

	if (M > N) {
		mnmax = M;
	} else {
		mnmax = N;
	}
	retcode = 0;
	P = 1;
	np = N + P;
	eps = macheps();

	AB = (double*) malloc(sizeof(double) * mnmax * np);
	q = (double*) malloc(sizeof(double) * N);
	C = (double*) malloc(sizeof(double) * N);
	V = (double*) malloc(sizeof(double) * N * N);

	for(i = 0; i < M;++i) {
		t = i * N;
		t2 = i * np;
		for(j = 0; j < N;++j) {
			AB[t2 + j] = Ai[t+j];
		}
	}

	for(i = 0; i < M;++i) {
		t2 = i * np;
		AB[t2+N] = bi[i];
	}

	minfit(AB,M,N,P,q);

	for(i = 0; i < N;++i) {
		t = i * N;
		t2 = i * np;
		for(j = 0; j < N;++j) {
			V[t+j] = AB[t2 + j];
		}
	}

	for(i = 0; i < N;++i) {
		t = i *N;
		for(j = 0; j < N;++j) {
			if (fabs(q[j]) > eps) {
				V[t+j] /= q[j];
			} else {
				V[t+j] = 0.0;
			}
		}
	}

	for(i = 0; i < N;++i) {
		t2 = i * np;
		C[i] = AB[t2+N];
	}
	mmult(V,C,xo,N,N,1);

	free(AB);
	free(q);
	free(V);
	free(C);
	return retcode;
}

int lls_svd(double *Ai,double *bi,int M,int N,double *xo) {
	int retcode,i,j,t;
	double *U,*V,*q,*C;
	double eps;

	U = (double*) malloc(sizeof(double) * M*N);
	V = (double*) malloc(sizeof(double) * N*N);
	q = (double*) malloc(sizeof(double) * N);
	C = (double*) malloc(sizeof(double) * N);

	retcode = 0;
	eps = macheps();

	retcode = svd(Ai,M,N,U,V,q);;
	if (retcode != 0) {
		return retcode;
	}

	mmult(bi,U,C,1,M,N);

	for(i = 0; i < N;++i) {
		t = i *N;
		for(j = 0; j < N;++j) {
			if (fabs(q[j]) > eps) {
				V[t+j] /= q[j];
			} else {
				V[t+j] = 0.0;
			}
		}
	}

	mmult(V,C,xo,N,N,1);

	free(U);
	free(V);
	free(q);

	return retcode;
}

int svd2( int m, int n, int withu, int withv, double eps, double tol, double a[4][3], double *q, double u[4][3], double v[3][3] ) {
	int i,j,k,l,l1,iter,retval;
	double c,f,g,h,s,x,y,z;
	double *e;

	e = (double *)calloc( n, sizeof(double) );
	retval = 0;

    // Copy 'a' to 'u'
	for ( i = 0; i < m; i++ ) {
		for ( j = 0; j < n; j++ )
			u[i][j] = a[i][j];
	}
    // Householder's reduction to bidiagonal form.
	g = x = 0.0;
	for (i = 0;i < n;i++) {
		e[i] = g;
		s = 0.0;
		l = i+1;
		for (j = i;j < m; j++ )
			s += ( u[j][i] * u[j][i] );
		if (s < tol)
			g = 0.0;
		else {
			f = u[i][i];
			g = (f < 0) ? sqrt(s) : -sqrt(s);
			h = f * g - s;
			u[i][i] = f - g;
			for ( j = l;j < n; j++ ) {
				s = 0.0;
				for ( k = i; k < m; k++ )
					s += ( u[k][i] * u[k][j] );
				f = s / h;
				for ( k=i; k<m; k++ )
					u[k][j] += (f * u[k][i]);
			} // end j
		} // end s
		q[i] = g;
		s = 0.0;
		for (j=l; j<n; j++)
			s += (u[i][j] * u[i][j]);
		if (s < tol)
			g = 0.0;
		else {
			f = u[i][i+1];
			g = (f < 0) ? sqrt(s) : -sqrt(s);
			h = f * g - s;
			u[i][i+1] = f - g;
			for (j=l;j<n;j++)
				e[j] = u[i][j]/h;
			for (j=l;j<m;j++) {
				s = 0.0;
				for (k=l;k<n;k++)
					s += (u[j][k] * u[i][k]);
				for (k=l;k<n;k++)
					u[j][k] += (s * e[k]);
			} // end j
		} // end s
		y = fabs(q[i]) + fabs(e[i]);
		if (y > x)
			x = y;
	} // end i



// accumulation of right-hand transformations
	if (withv) {
		for( i = n - 1; i >= 0; i--) {
		   if ( i < n - 1 ) {                   //============================gg
			  if ( g != 0.0)  {
				 h = u[ i ][ i + 1 ]  *  g;     //==============in other implementation: h[i][l]
				 for( j = l; j < n; j++ )
					v[ j ][ i ] = u[ i ][ j ] / h;
				 for( j = l; j < n; j++ ) {
					s = 0.0;
					for (k = l; k < n; k++ )
						s += ( u[ i ][ k ] * v[ k ][ j ] );
					for (k = l; k < n; k++ )
						v[ k ][ j ] += ( s * v[ k ][ i ] );

				 } // end j
			  } // end g
			  for ( j = l; j < n; j++ )
				  v[ i ][ j ] = v[ j ][ i ] = 0.0;
			}                                   // end if(i < n-1) ===========gg
		    v[ i ][ i ] = 1.0;
			g = e[ i ];
			l = i;
		} // end i

	} // end withv, parens added for clarity

	if ( withu ) {
		for ( i = n - 1; i >= 0; i-- ) {
			l = i + 1;
			g = q[ i ];
			if ( i < n - 1 )                    //============================gg
			   for( j = l; j < n; j++ )         // upper limit was 'n' (--> 'm' in original===gg)
				  u[ i ][ j ] = 0.0;
			if ( g != 0.0 ) {
			   if ( i != n - 1 ) {              //============================gg
			      h = u[ i ][ i ] * g;
			      for( j = l; j < n; j++ ) {    // upper limit was 'n' (--> 'm' in original===gg)
				     s = 0.0;
				     for( k = l; k < m; k++ )
						s += ( u[ k ][ i ] * u[ k ][ j ] );
				     f = s / h;
				     for( k = i; k < m; k++ )
					    u[ k ][ j ] += ( f * u[ k ][ i ] );
				   }  // end j
				}     // end if( i != n-1 )      =============================gg
				for ( j = i; j < m; j++ )
					u[ j ][ i ] /= g;
			  } // end if (g)
			  else {
				for ( j = i; j < m; j++)
					u[ j ][ i ] = 0.0;
			}
			u[ i ][ i ] += 1.0;
		} // end i
	} // end withu, parens added for clarity

// diagonalization of the bidiagonal form
	for(i = 0; i < m;++i) {
		for(j = 0; j < n;++j) {
			printf("%g ",u[i][j]);
		}
		printf("\n");
	}

	eps *= x;
	for (k = n - 1; k >= 0; k--) {
		iter = 0;
test_f_splitting:
		for (l = k ; l >= 0;l-- ) {
			if (fabs(e[l]) <= eps) goto test_f_convergence;
			if (fabs(q[l-1]) <= eps) goto cancellation;
		} // end l

// cancellation of e[l] if l > 0
cancellation:
		c = 0.0;
		s = 1.0;
		l1 = l - 1;
		for (i = l; i <= k; i++) {
			f = s * e[ i ];
			e[ i ] *= c;
			if ( fabs(f) <= eps ) goto test_f_convergence;
			g = q[ i ];
			h = q[ i ] = sqrt( f*f + g*g );
			c = g / h;
			s = -f / h;
			if ( withu ) {
			   for( j = 0; j < m; j++ ) {
					y = u[ j ][ l1 ];
					z = u[ j ][ i ];
					u[ j ][ l1 ] = y * c + z * s;
					u[ j ][ i ] = -y * s + z * c;
				} // end j
			} // end withu, parens added for clarity
		} // end i
test_f_convergence:
		z = q[ k ];
		if (l == k) goto convergence;

// shift from bottom 2x2 minor
		iter++;
		if (iter > 30) {
			retval = k;
			break;
		}
		x = q[l];
		y = q[k-1];
		g = e[k-1];
		h = e[k];
		f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2*h*y);
		g = sqrt(f*f + 1.0);
		f = ((x-z)*(x+z) + h*(y/((f<0)?(f-g):(f+g))-h))/x;
// next QR transformation
		c = s = 1.0;
		for (i=l+1;i<=k;i++) {
			g = e[i];
			y = q[i];
			h = s * g;
			g *= c;
			e[i-1] = z = sqrt(f*f+h*h);
			c = f / z;
			s = h / z;
			f = x * c + g * s;
			g = -x * s + g * c;
			h = y * s;
			y *= c;
			if (withv) {
				for (j=0;j<n;j++) {
					x = v[j][i-1];
					z = v[j][i];
					v[j][i-1] = x * c + z * s;
					v[j][i] = -x * s + z * c;
				} // end j
			} // end withv, parens added for clarity
			q[i-1] = z = sqrt(f*f + h*h);
			c = f/z;
			s = h/z;
			f = c * g + s * y;
			x = -s * g + c * y;
			if (withu) {
				for (j=0;j<m;j++) {
					y = u[j][i-1];
					z = u[j][i];
					u[j][i-1] = y * c + z * s;
					u[j][i] = -y * s + z * c;
				} // end j
			} // end withu, parens added for clarity
		} // end i
		e[l] = 0.0;
		e[k] = f;
		q[k] = x;
		goto test_f_splitting;
convergence:
		if (z < 0.0) {
// q[k] is made non-negative
			q[k] = - z;
			if (withv) {
				for (j=0;j<n;j++)
					v[j][k] = -v[j][k];
			} // end withv, parentheses added for clarity
		} // end z
	} // end k

	free(e);
	return retval;
}

void svd_csa(double A[4][3], int m,int n, double* w, double V[3][3])
{
    double* rv1;
    int i, j, k, l = -1;
    double tst1, c, f, g, h, s, scale;


    rv1 = malloc(n * sizeof(double));

    /*
     * householder reduction to bidiagonal form
     */

    g = 0.0;
    scale = 0.0;
    tst1 = 0.0;
    for (i = 0; i < n; i++) {



        l = i + 1;
        rv1[i] = scale * g;
        g = 0.0;
        s = 0.0;
        scale = 0.0;
        if (i < m) {
            for (k = i; k < m; k++)
                scale += fabs(A[k][i]);
            if (scale != 0.0) {
                for (k = i; k < m; k++) {
                    A[k][i] /= scale;
                    s += A[k][i] * A[k][i];
                }
                f = A[i][i];
                g = -copysign(sqrt(s), f);
                h = f * g - s;
                A[i][i] = f - g;
                if (i < n - 1) {        /* no test in NR */
                    for (j = l; j < n; j++) {
                        s = 0.0;
                        for (k = i; k < m; k++)
                            s += A[k][i] * A[k][j];
                        f = s / h;
                        for (k = i; k < m; k++)
                            A[k][j] += f * A[k][i];
                    }
                }
                for (k = i; k < m; k++)
                    A[k][i] *= scale;
            }
        }
        w[i] = scale * g;
        g = 0.0;
        s = 0.0;
        scale = 0.0;
        if (i < m && i < n - 1) {
            for (k = l; k < n; k++)
                scale += fabs(A[i][k]);
            if (scale != 0.0) {
                for (k = l; k < n; k++) {
                    A[i][k] /= scale;
                    s += A[i][k] * A[i][k];
                }
                f = A[i][l];
                g = -copysign(sqrt(s), f);
                h = f * g - s;
                A[i][l] = f - g;
                for (k = l; k < n; k++)
                    rv1[k] = A[i][k] / h;
                for (j = l; j < m; j++) {
                    s = 0.0;
                    for (k = l; k < n; k++)
                        s += A[j][k] * A[i][k];
                    for (k = l; k < n; k++)
                        A[j][k] += s * rv1[k];
                }
                for (k = l; k < n; k++)
                    A[i][k] *= scale;
            }
        }
        {
            double tmp = fabs(w[i]) + fabs(rv1[i]);

            tst1 = (tst1 > tmp) ? tst1 : tmp;
        }
    }

	for(i = 0; i < m;++i) {
		for(j = 0; j < n;++j) {
			printf("%g ",A[i][j]);
		}
		printf("\n");
	}

    /*
     * accumulation of right-hand transformations
     */

    for (i = n - 1; i >= 0; i--) {


        if (i < n - 1) {        /* no test in NR */
            if (g != 0.0) {
                for (j = l; j < n; j++)
                    /*
                     * double division avoids possible underflow
                     */
                    V[j][i] = (A[i][j] / A[i][l]) / g;
                for (j = l; j < n; j++) {
                    s = 0.0;
                    for (k = l; k < n; k++)
                        s += A[i][k] * V[k][j];
                    for (k = l; k < n; k++)
                        V[k][j] += s * V[k][i];
                }
            }
            for (j = l; j < n; j++) {
                V[i][j] = 0.0;
                V[j][i] = 0.0;
            }
        }
        V[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /*
     * accumulation of left-hand transformations
     */

    for (i = (m < n) ? m - 1 : n - 1; i >= 0; i--) {



        l = i + 1;
        g = w[i];
        if (i != n - 1)
            for (j = l; j < n; j++)
                A[i][j] = 0.0;
        if (g != 0.0) {
            for (j = l; j < n; j++) {
                s = 0.0;
                for (k = l; k < m; k++)
                    s += A[k][i] * A[k][j];
                /*
                 * double division avoids possible underflow
                 */
                f = (s / A[i][i]) / g;
                for (k = i; k < m; k++)
                    A[k][j] += f * A[k][i];
            }
            for (j = i; j < m; j++)
                A[j][i] /= g;
        } else
            for (j = i; j < m; j++)
                A[j][i] = 0.0;
        A[i][i] += 1.0;
    }

    /*
     * diagonalization of the bidiagonal form
     */

    for (k = n - 1; k >= 0; k--) {
        int k1 = k - 1;
        int its = 0;


        while (1) {
            int docancellation = 1;
            double x, y, z;
            int l1 = -1;

            its++;
            if (its > SVDMAXITER)
                exit(1);

            for (l = k; l >= 0; l--) {  /* test for splitting */
                double tst2 = fabs(rv1[l]) + tst1;

                if (tst2 == tst1) {
                    docancellation = 0;
                    break;
                }
                l1 = l - 1;
                /*
                 * rv1(1) is always zero, so there is no exit through the
                 * bottom of the loop
                 */
                tst2 = fabs(w[l - 1]) + tst1;
                if (tst2 == tst1)
                    break;
            }
            /*
             * cancellation of rv1[l] if l > 1
             */
            if (docancellation) {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) {
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if ((fabs(f) + tst1) == tst1)
                        break;
                    g = w[i];
                    h = hypot(f, g);
                    w[i] = h;
                    c = g / h;
                    s = -f / h;
                    for (j = 0; j < m; j++) {
                        double y = A[j][l1];
                        double z = A[j][i];

                        A[j][l1] = y * c + z * s;
                        A[j][i] = z * c - y * s;
                    }
                }
            }
            /*
             * test for convergence
             */
            z = w[k];
            if (l != k) {
                int i1;

                /*
                 * shift from bottom 2 by 2 minor
                 */
                x = w[l];
                y = w[k1];
                g = rv1[k1];
                h = rv1[k];
                f = 0.5 * (((g + z) / h) * ((g - z) / y) + y / h - h / y);
                g = hypot(f, 1.0);
                f = x - (z / x) * z + (h / x) * (y / (f + copysign(g, f)) - h);
                /*
                 * next qr transformation
                 */
                c = 1.0;
                s = 1.0;
                for (i1 = l; i1 < k; i1++) {
                    i = i1 + 1;
                    g = rv1[i];
                    y = w[i];
                    h = s * g;
                    g = c * g;
                    z = hypot(f, h);
                    rv1[i1] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = g * c - x * s;
                    h = y * s;
                    y *= c;
                    for (j = 0; j < n; j++) {
                        x = V[j][i1];
                        z = V[j][i];
                        V[j][i1] = x * c + z * s;
                        V[j][i] = z * c - x * s;
                    }
                    z = hypot(f, h);
                    w[i1] = z;
                    /*
                     * rotation can be arbitrary if z = 0
                     */
                    if (z != 0.0) {
                        c = f / z;
                        s = h / z;
                    }
                    f = c * g + s * y;
                    x = c * y - s * g;
                    for (j = 0; j < m; j++) {
                        y = A[j][i1];
                        z = A[j][i];
                        A[j][i1] = y * c + z * s;
                        A[j][i] = z * c - y * s;
                    }
                }
                rv1[l] = 0.0;
                rv1[k] = f;
                w[k] = x;
            } else {
                /*
                 * w[k] is made non-negative
                 */
                if (z < 0.0) {
                    w[k] = -z;
                    for (j = 0; j < n; j++)
                        V[j][k] = -V[j][k];
                }
                break;
            }
        }
    }


    free(rv1);
}




