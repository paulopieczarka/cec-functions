#include "cec19_func.hpp"

Cec19::Cec19 () :
	ini_flag(0),
	M(nullptr),
	OShift(nullptr),
	z(nullptr),
	y(nullptr),
	x_bound(nullptr)
	{}

Cec19::~Cec19 () {
	this->free_func();
}

void Cec19::test_func (double *x, double *f, int nx, int mx,int func_num) {
	int cf_num = 10, i, j;

	if (ini_flag == 1) {
		if ((n_flag != nx) || (func_flag != func_num)) {
			ini_flag=0;
		}
	}

	if (ini_flag == 0) {
		FILE *fpt;
		char FileName[256];

		this->free_func();

		y = (double *)malloc(sizeof(double)  *  nx);
		z = (double *)malloc(sizeof(double)  *  nx);
		x_bound = (double *)malloc(sizeof(double)  *  nx);

		for (i=0; i<nx; i++) {
			x_bound[i] = 100.0;
		}

		if (!(nx==2 || nx == 10 || nx == 9 || nx == 16 || nx == 18)) {
			std::cout << "Error: Test functions are only defined for D=10, 9, 16, 18" << std::endl;
			std::cout << "F1 is defined on D=9 \n F2 is defined on D=16" << std::endl;
			std::cout << "F3 is defined on D=18 \n F4-F10 are defined on D=10." << std::endl;
		}

		/* Load Matrix M*/
		if (func_num > 3) {
			sprintf(FileName, "input_data/M_%d_D%d.txt", func_num,nx);
			fpt = fopen(FileName,"r");

			if (fpt == NULL) {
				std::cout << "Error: Cannot open input file for reading " << std::endl;
			}

			M = (double*)malloc(nx * nx * sizeof(double));
			if (M == NULL) {
				std::cout << "Error: there is insufficient memory available!" << std::endl;
			}

			for (i = 0; i < (nx * nx); i++) {
				fscanf(fpt, "%lf", &M[i]);
			}

			fclose(fpt);
		}

		/* Load shift_data */
		if (func_num > 3) {
			sprintf(FileName, "input_data/shift_data_%d.txt", func_num);
			fpt = fopen(FileName,"r");
			if (fpt==NULL) {
				std::cout << "Error: Cannot open input file for reading" << std::endl;
			}

			OShift = (double *)malloc(nx * sizeof(double));
			if (OShift == NULL) {
				std::cout << "Error: there is insufficient memory available!" << std::endl;
			}

			for(i = 0; i < nx; i++) {
				fscanf(fpt, "%lf", &OShift[i]);
			}

			fclose(fpt);
		}

		n_flag = nx;
		func_flag = func_num;
		ini_flag = 1;
	}

	for (i = 0; i < mx; i++) {
		switch(func_num) {
			case 1:
				Chebyshev_func(&x[i*nx], nx, &f[i]);
				f[i] += 1.0;
				break;

			case 2:
				Hilbert_func(&x[i*nx], nx, &f[i]);
				f[i] += 1.0;
				break;

			case 3:
				Lennard_Jones_func(&x[i*nx], nx, &f[i]);
				f[i] += 1.0;
				break;

			case 4:
				Rastrigin_func(&x[i*nx], &f[i], nx, OShift, M, 1, 1);
				f[i] += 1.0;
				break;

			case 5:
				Griewank_func(&x[i*nx], &f[i], nx, OShift, M, 1, 1);
				f[i] += 1.0;
				break;

			case 6:
				Weierstrass_func(&x[i*nx], &f[i], nx, OShift, M, 1, 1);
				f[i] += 1.0;
				break;

			case 7:
				Schwefel_func(&x[i*nx], &f[i], nx, OShift, M, 1, 1);
				f[i] += 1.0;
				break;

			case 8:
				EScaffer6_func(&x[i*nx], &f[i], nx, OShift, M, 1, 1);
				f[i] += 1.0;
				break;

			case 9:
				Happycat_func(&x[i*nx], &f[i], nx, OShift, M, 1, 1);
				f[i] += 1.0;
				break;

			case 10:
				Ackley_func(&x[i*nx], &f[i], nx, OShift, M, 1, 1);
				f[i] += 1.0;
				break;

			default:
				std::cout << "Error: There are only 30 test functions in this test suite!" << std::endl;
				f[i] = 0.0;
				break;
		}
	}

}
/* Schwefel's 1.2 */
void Cec19::Schaffer_F7_func (double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag) {
  int i;
	double tmp;
  f[0] = 0.0;

	// shift and rotate
	sr_func (x, this->z, nx, Os, Mr, 1.0, s_flag, r_flag);

	for (i=0; i < nx-1; i++) {
		this->z[i] = pow(this->y[i] *y[i] +y[i+1] * y[i+1], 0.5);
		tmp=sin(50.0 * pow(this->z[i], 0.2));
		f[0] += pow(this->z[i], 0.5) + pow(this->z[i], 0.5) * tmp * tmp;
	}

	f[0] = f[0] * f[0] / (nx-1) / (nx-1);
}

/* Griewank's */
void Cec19::Griewank_func (double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag) {
  int i;
  double s, p;
  s = 0.0;
  p = 1.0;

	// shift and rotate
	sr_func (x, this->z, nx, Os, Mr, 600.0/100.0, s_flag, r_flag);

	for (i=0; i<nx; i++) {
		s += this->z[i] * this->z[i];
		p *= cos(this->z[i] / sqrt(1.0+i));
	}

	f[0] = 1.0 + s / 4000.0 - p;
}

/* Ackley's */
void Cec19::Ackley_func (double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag) {
  int i;
  double sum1, sum2;
  sum1 = 0.0;
  sum2 = 0.0;

	// shift and rotate
	sr_func (x, this->z, nx, Os, Mr, 1.0, s_flag, r_flag);

	for (i=0; i<nx; i++) {
		sum1 += this->z[i] * this->z[i];
		sum2 += cos(2.0 * PI * this->z[i]);
	}

	sum1 = -0.2 * sqrt(sum1 / nx);
	sum2 /= nx;

	f[0] =  E - 20.0*exp(sum1) - exp(sum2) +20.0;
}

/* Weierstrass's */
void Cec19::Weierstrass_func (double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag) {
    int i, j, k_max;
    double sum, sum2, a, b;

    a = 0.5;
    b = 3.0;
    k_max = 20;
    f[0] = 0.0;

	// shift and rotate
	sr_func (x, this->z, nx, Os, Mr, 0.5/100.0, s_flag, r_flag);

	for (i = 0; i < nx; i++) {
		sum = 0.0;
		sum2 = 0.0;
		for (j=0; j <= k_max; j++) {
			sum += pow(a, j) * cos(2.0 * PI * pow(b, j) * (this->z[i] + 0.5));
			sum2 += pow(a, j) * cos(2.0 * PI * pow(b, j) * 0.5);
		}

		f[0] += sum;
	}

	f[0] -= nx * sum2;
}

/* Rastrigin's */
void Cec19::Rastrigin_func (double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag) {
  int i;
	f[0] = 0.0;

	// shift and rotate
	sr_func (x, this->z, nx, Os, Mr, 5.12/100.0, s_flag, r_flag);

	for (i=0; i<nx; i++) {
		f[0] += (this->z[i] * this->z[i] - 10.0 * cos(2.0 * PI * this->z[i]) + 10.0);
	}
}

/* Noncontinuous Rastrigin's */
void Cec19::Step_rastrigin_func (double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag) {
  int i;
	f[0] = 0.0;

	for (i = 0; i < nx; i++) {
		if (fabs(this->y[i] - Os[i]) > 0.5) {
			this->y[i] = Os[i] + floor(2 * (this->y[i] - Os[i]) + 0.5) / 2;
		}
	}

	// shift and rotate
	sr_func (x, this->z, nx, Os, Mr, 5.12 / 100.0, s_flag, r_flag);

	for (i=0; i < nx; i++) {
		f[0] += (this->z[i] * this->z[i] - 10.0 * cos(2.0 * PI * this->z[i]) + 10.0);
	}
}

/* Schwefel's */
void Cec19::Schwefel_func (double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag) {
  int i;
	double tmp;
	f[0] = 0.0;

	// shift and rotate
	sr_func (x, this->z, nx, Os, Mr, 1000.0/100.0, s_flag, r_flag);

	for (i=0; i < nx; i++) {
		this->z[i] += 4.209687462275036e+002;

		if (this->z[i] > 500) {
			f[0] -= (500.0 - fmod(this->z[i], 500)) * sin(pow(500.0 - fmod(this->z[i], 500), 0.5));
			tmp=(this->z[i]-500.0)/100;
			f[0]+= tmp*tmp/nx;
		} else if (this->z[i] <= 500) {
			f[0] -= (-500.0 + fmod(fabs(this->z[i]), 500)) * sin(pow(500.0 - fmod(fabs(this->z[i]), 500), 0.5));
			tmp = (this->z[i] + 500.0) / 100;
			f[0] += tmp * tmp / nx;
		} else {
			f[0] -= this->z[i]*sin(pow(fabs(this->z[i]),0.5));
		}
	}

	f[0] +=4.189828872724338e+002*nx;
}

/* Expanded Scaffer's F6 */
void Cec19::EScaffer6_func (double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag) {
  int i;
  double temp1, temp2;

	// shift and rotate
	sr_func (x, this->z, nx, Os, Mr, 1.0, s_flag, r_flag);

  f[0] = 0.0;
  for (i=0; i<nx-1; i++) {
    temp1 = sin(sqrt(this->z[i] * this->z[i] + this->z[i+1] * this->z[i+1]));
		temp1 = temp1 * temp1;
    temp2 = 1.0 + 0.001 * (this->z[i] * this->z[i] + this->z[i+1] * this->z[i+1]);
    f[0] += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
  }

	temp1 = sin(sqrt(this->z[nx-1] * this->z[nx-1] + this->z[0] * this->z[0]));
	temp1 = temp1 * temp1;
  temp2 = 1.0 + 0.001 * (this->z[nx-1] * this->z[nx-1] + this->z[0] * this->z[0]);
  f[0] += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
}

/* HappyCat, provdided by Hans-Georg Beyer (HGB) */
void Cec19::Happycat_func (double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag) {
	int i;
	double alpha, r2, sum_z;
	alpha = 1.0 / 8.0;

	// shift and rotate
	sr_func (x, this->z, nx, Os, Mr, 5.0/100.0, s_flag, r_flag);

	r2 = 0.0;
	sum_z = 0.0;

	for (i=0; i<nx; i++) {
		this->z[i] = this->z[i] - 1.0; //shift to orgin
    r2 += this->z[i] * this->z[i];
		sum_z += this->z[i];
  }

	f[0] = pow(fabs(r2-nx), 2 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5;
}

void Cec19::shift_func (double *x, double *xshift, int nx, double *Os) {
  for (int i = 0; i < nx; i++) {
    xshift[i] = x[i] - Os[i];
  }
}

void Cec19::rotate_func (double *x, double *xrot, int nx, double *Mr) {
  for (int i = 0; i < nx; i++) {
    xrot[i] = 0;
		for (int j = 0; j < nx; j++) {
			xrot[i] = xrot[i] + x[j] * Mr[i * nx + j];
		}
  }
}

void Cec19::sr_func (double *x, double *sr_x, int nx, double *Os, double *Mr, double sh_rate, int s_flag, int r_flag) {
	int i;
	if (s_flag == 1) {
		if (r_flag == 1) {
			this->shift_func(x, y, nx, Os);

			//shrink to the orginal search range
			for (i=0; i<nx; i++) {
				y[i] = y[i] * sh_rate;
			}

			this->rotate_func(y, sr_x, nx, Mr);
		} else {
			this->shift_func(x, sr_x, nx, Os);

			//shrink to the orginal search range
			for (i=0; i<nx; i++) {
				sr_x[i] = sr_x[i] * sh_rate;
			}
		}
	} else {
		if (r_flag == 1) {
			//shrink to the orginal search range
			for (i = 0; i < nx; i++) {
				y[i] = x[i] * sh_rate;
			}

			this->rotate_func(y, sr_x, nx, Mr);
		} else {
			//shrink to the orginal search range
			for (i = 0; i < nx; i++) {
				sr_x[i] = x[i] * sh_rate;
			}
		}
	}
}

void Cec19::asy_func (double *x, double *xasy, int nx, double beta) {
  for (int i = 0; i < nx; i++) {
		if (x[i] > 0) {
      xasy[i] = pow(x[i], 1.0 + beta * i / (nx - 1) * pow(x[i], 0.5));
		}
  }
}

void Cec19::osz_func (double *x, double *xosz, int nx) {
	int i, sx;
	double c1, c2, xx;

	for (i = 0; i < nx; i++) {
		if (i == 0 || i == nx-1) {
			if (x[i] != 0) {
				xx = log(fabs(x[i]));
			}

			if (x[i] > 0) {
				c1 = 10;
				c2 = 7.9;
			} else {
				c1 = 5.5;
				c2 = 3.1;
			}

			if (x[i] > 0) {
				sx = 1;
			} else if (x[i] == 0) {
				sx = 0;
			} else {
				sx = -1;
			}

			xosz[i] = sx * exp(xx + 0.049 * (sin(c1 * xx) + sin(c2 * xx)));
		}
		else {
			xosz[i]=x[i];
		}
  }
}



/**
	find the atomic configuration with minimum energy

	valid for any dimension, D=3*k, k=2,3,4,...,25.   k is the number of atoms in 3-D space
	constraints: unconstrained
	type: multi-modal with one global minimum; non-separable
	initial upper bound = 4, initial lower bound = -4
	value-to-reach = minima[k-2]+.0001
	f(x*) = minima[k-2]; see array of minima below; additional minima available at the
	Cambridge cluster database: http://www-wales.ch.cam.ac.uk/~jon/structures/LJ/tables.150.html
*/
void Cec19::Lennard_Jones_func (double *x, int D, double *f) {
	f[0] = 0;
	int i, j, k, a, b;
	long double xd, yd, zd, ed, ud, sum = 0;

	static double minima[] = {
		-1.,-3.,-6., -9.103852, -12.712062, -16.505384, -19.821489, -24.113360,
		-28.422532, -32.765970, -37.967600, -44.326801, -47.845157, -52.322627, -56.815742, -61.317995,
		-66.530949, -72.659782, -77.1777043, -81.684571, -86.809782, -02.844472, -97.348815, -102.372663
	};

	k = D / 3;

	// default if k<2
	if (k < 2) {
		k = 2;
		D = 6;
	}

	for (i = 0; i < k - 1; i++) {
		for (j = i + 1; j < k; j++) {
			a = 3 * i;
			b = 3 * j;
			xd = x[a] - x[b];
			yd = x[a + 1] - x[b + 1];
			zd = x[a + 2] - x[b + 2];
			ed = xd*xd + yd*yd + zd*zd;
			ud = ed*ed*ed;
			if (ud > 1.0e-10) sum += (1.0 / ud - 2.0) / ud;
			else sum += 1.0e20;
		}
	}

	f[0] += sum;
	f[0] += 12.7120622568;
}

/**
	find the inverse of the (ill-conditioned) Hilbert matrix

	valid for any dimension, n=k*k, k=2,3,4,...
	constraints: unconstrained
	type: multi-modal with one global minimum; non-separable
	initial upper bound = 2^n, initial lower bound = -2^n
	value-to-reach = f(x*)+1.0e-8
	f(x*) = 0.0; x*={{9,-36,30},{-36,192,-180},{30,-180,180}} (n=9)
	x*={{16,-120,240,-140},{-120,1200,-2700,1680},{240,-2700,6480,4200},{-140,1680,-4200,2800}} (n=16)
*/
void Cec19::Hilbert_func (double *x, int D, double *f) {
	f[0] = 0;
	int i, j, k, b;

	long double sum = 0;

	// Increase matrix size if D > 100
	static long double hilbert[10][10], y[10][10];

	b = (int)sqrt((double)D);

	for (i = 0; i < b; i++) {
		for (j = 0; j < b; j++) {
			// Create a static Hilbert matrix
			hilbert[i][j] = 1. / (double)(i + j + 1);
		}
	}

	for (j = 0; j < b; j++) {
		for (k = 0; k < b; k++) {
			y[j][k] = 0;
			for (i = 0; i < b; i++) {
				// Compute matrix product H*x
				y[j][k] += hilbert[j][i] * x[k + b * i];
			}
		}
	}

	for (i = 0; i < b; i++) {
		for (j = 0; j < b; j++) {
			// Sum absolute value of deviations
			if (i == j) {
				sum += fabs(y[i][j] - 1);
			} else {
				sum += fabs(y[i][j]);
			}
		}
	}

	f[0] += sum;
}

/**
	Storn's Tchebychev - a 2nd ICEO function - generalized version

	Valid for any D>2
	constraints: unconstrained
	type: multi-modal with one global minimum; non-separable
	initial upper bound = 2^D, initial lower bound = -D^n
	value-to-reach = f(x*)+1.0e-8
	f(x*)=0.0; x*=(128,0,-256,0,160,0,-32,0,1) (n=9)
	x*=(32768,0,-131072,0,212992,0,-180224,0,84480,0,-21504,0,2688,0,-128,0,1) (n=17)
*/
void Cec19::Chebyshev_func (double *x, int D, double *f) {
	f[0] = 0.0;
	int i, j;
	static int sample;
	long double a = 1.0, b = 1.2, px, y = -1, sum = 0;
	static long double dx, dy;

	for (j = 0; j < D - 2; j++) {
		dx = 2.4 * b - a;
		a = b;
		b = dx;
	}

	sample = 32 * D;
	dy = 2.0 / (long double)sample;

	for (i = 0; i <= sample; i++) {
		px = x[0];
		for (j = 1; j < D; j++) {
			px = y*px + x[j];
		}

		if (px < -1 || px > 1) {
			sum += (1. - fabs(px))*(1. - fabs(px));
		}

		y += dy;
	}

	for (i = -1; i <= 1; i += 2) {
		px = x[0];
		for (j = 1; j < D; j++) {
			px = 1.2 * px + x[j];
		}

		if (px < dx) {
			sum += px * px;
		}
	}

	f[0] += sum;
}

void Cec19::free_func () {
	std::free(this->M);
	std::free(this->OShift);
	std::free(this->y);
	std::free(this->z);
	std::free(this->x_bound);
	z = nullptr;
  M = nullptr;
  OShift = nullptr;
  x_bound = nullptr;
  ini_flag = 0;
}
