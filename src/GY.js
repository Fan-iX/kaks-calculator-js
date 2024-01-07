/************************************************************
* Filename: GY.js
* Abstract: Definition of GY class.

* Version: js-3.1
* Author: Fanix
* Date: January.4, 2024

* Re-implementation of KaKs_Calculator/GY94

  References:
  Goldman N, Yang Z  (1994)  A codon-based model of nucleotide
  substitution for protein-coding DNA sequences. Mol. Biol. 
  Evol. 11:725-736.
*************************************************************/

class GY extends Base {
	/*
	JC,    F81:   rTC==rAG =rTA==rCG==rTG==rCA
	K2P,   HKY:   rTC==rAG!=rTA==rCG==rTG==rCA
	TNEF,  TN:    rTC!=rAG!=rTA==rCG==rTG==rCA
	K3P,   K3PUF: rTC==rAG!=rTA==rCG!=rTG==rCA
	TIMEF, TIM:   rTC!=rAG!=rTA==rCG!=rTG==rCA
	TVMEF, TVM:   rTC==rAG!=rTA!=rCG!=rTG!=rCA
	SYM,   GTR:   rTC!=rAG!=rTA!=rCG!=rTG!=rCA
	*/
	constructor(genetic_code = 1, model = "HKY") {
		super(genetic_code)
		this.method = "GY-" + model;
		this.model = model;
		this.nkappa = {
			"JC": 0, "F81": 0,
			"K2P": 1, "HKY": 1,
			"TNEF": 2, "TN": 2,
			"K3P": 2, "K3PUF": 2,
			"TIMEF": 3, "TIM": 3,
			"TVMEF": 4, "TVM": 4,
			"SYM": 5, "GTR": 5
		}[this.model];
		if (this.nkappa == null) throw "Invalid model";

		this.FROM61 = VectorUtils.seq(this.translation.length).filter(i => this.translation[i] != "!");
		this.FROM64 = [];
		let k = 0;
		this.translation.split("").forEach((x, i) => { this.FROM64[i] = x != "!" ? k : -1; if (x != "!") k++ })
	}

	/* x[i]=x0[i] + a*p[i] */
	fun_ls(a, x0, p) {
		let x = x0.map((_, i) => x0[i] + a * p[i])
		return this.lfun2dSdN(x);
	}

	/**
	* Count codon frequence
	* @param {string} seq encoded DNA sequence
	* @returns {number[]} `pi64` codon frequence
	*/
	getFreqency(seq) {
		let pi = new Array(64).fill(0);
		let arr = seq.match(/.{3}/g);
		if (
			new Set(arr.filter(c => this.getAminoAcid(c) != "!")).size == this.translation.split("").filter(x => x != "!").length
		) {
			for (let c of arr) pi[this.encodeCodon(c)]++
		} else {
			//Get A,C,G,T frequency at three positions
			let f12pos = new Array(12).fill(0);
			seq.split("").forEach((x, i) => f12pos[(i % 3) * 4 + +x]++)
			f12pos = f12pos.map(x => x / (seq.length / 3))

			// Get 64 amino acid probability
			for (let i of [0, 1, 2, 3])
				for (let j of [0, 1, 2, 3])
					for (let k of [0, 1, 2, 3])
						if (this.translation[i * 16 + j * 4 + k] != '!')
							pi[i * 16 + j * 4 + k] = f12pos[i] * f12pos[j + 4] * f12pos[k + 8];
		}
		//Scale the sum of pi[] to 1
		let scale = VectorUtils.sum(pi)
		pi = pi.map(x => x / scale)

		if (Math.abs(1 - VectorUtils.sum(pi)) > 1e-6)
			console.warn("Warning: error in get codon freqency.");
		return pi
	}

	gradientB(x, f0, xmark, sizep) {
		/* f0=fun(x) is always provided.
		   xmark=0: central; 1: upper; -1: down
		*/
		let eh0 = 1e-6, eh;  /* eh0=1e-6 || 1e-7 */
		let g = []
		for (let i = 0; i < x.length; i++) {
			eh = eh0 * (Math.abs(x[i]) + 1);
			if (xmark[i] == 0 && sizep < 1) {    //central 
				eh = Math.pow(eh, .67);
				let t1 = [...x], t2 = [...x]
				t1[i] -= eh; t2[i] += eh;
				g[i] = (this.lfun2dSdN(t2) - this.lfun2dSdN(t1)) / (eh * 2);
			}
			else {//forward or backward
				let t2 = [...x]
				if (xmark[i]) eh *= -xmark[i];
				t2[i] += eh;
				g[i] = (this.lfun2dSdN(t2) - f0) / eh;
			}
		}
		return g
	}

	LineSearch2(f, x0, p, step, limit, e) {
		/* linear search using quadratic interpolation 
		   from x0[] in the direction of p[],
						x = x0 + a*p        a ~(0,limit)
		   returns (a).    *f: f(x0) for input and f(x) for output
		
		   x0[n] x[n] p[n] space[n]
		
		   adapted from Wolfe M. A.  1978.  Numerical methods for unconstrained
		   optimization: An introduction.  Van Nostrand Reinhold Company, New York.
		   pp. 62-73.
		   step is used to find the bracket and is increased or reduced as necessary, 
		   and is not terribly important.
		*/
		let maxround = 10;
		let factor = 4, small = 1e-6, smallgapa = 0.2;
		let a3, a4, a5, a6, f3, f4, f5, f6;

		let a0 = 0, a1 = 0, f1 = f, f0 = f,
			a2 = a0 + step,
			f2 = this.fun_ls(a2, x0, p);
		if (f2 > f1) {
			while (f2 > f1) {
				step /= factor;
				if (step < small) return 0;
				a3 = a2; f3 = f2;
				a2 = a0 + step; f2 = this.fun_ls(a2, x0, p);
			}
		}
		else {
			for (; ;) {
				step *= factor;
				if (step > limit) step = limit;
				a3 = a0 + step; f3 = this.fun_ls(a3, x0, p);

				//obtain a bracket
				if (f3 >= f2) {	//a1<a2<a3 and f1>f2<f3
					break;
				}
				if (step >= limit) return a3;
				a1 = a2; f1 = f2; a2 = a3; f2 = f3;
			}
		}

		/* iteration by quadratic interpolation, fig 2.2.9-10 (pp 71-71) */
		for (let i = 0; i < maxround; i++) {	//a4 is the minimum from the parabola over (a1,a2,a3)

			//2.2.4
			a4 = (a2 - a3) * f1 + (a3 - a1) * f2 + (a1 - a2) * f3;
			if (Math.abs(a4) > 1e-100) {
				//2.2.3
				a4 = ((a2 * a2 - a3 * a3) * f1 + (a3 * a3 - a1 * a1) * f2 + (a1 * a1 - a2 * a2) * f3) / (2 * a4);
			}
			if (a4 > a3 || a4 < a1) //out of range: whether a1<a4<a3
				a4 = (a1 + a2) / 2;

			f4 = this.fun_ls(a4, x0, p);

			if (Math.abs(f2 - f4) < e * (1 + Math.abs(f2)))
				break;


			if (a4 <= a2) {    /* fig 2.2.10 */
				if (a2 - a4 > smallgapa * (a2 - a1)) {
					if (f4 <= f2) { a3 = a2; a2 = a4; f3 = f2; f2 = f4; }
					else { a1 = a4; f1 = f4; }
				}
				else {
					if (f4 > f2) {
						a5 = (a2 + a3) / 2; f5 = this.fun_ls(a5, x0, p);
						if (f5 > f2) { a1 = a4; a3 = a5; f1 = f4; f3 = f5; }
						else { a1 = a2; a2 = a5; f1 = f2; f2 = f5; }
					}
					else {
						a5 = (a1 + a4) / 2; f5 = this.fun_ls(a5, x0, p);
						if (f5 >= f4) { a3 = a2; a2 = a4; a1 = a5; f3 = f2; f2 = f4; f1 = f5; }
						else {
							a6 = (a1 + a5) / 2; f6 = this.fun_ls(a6, x0, p);
							if (f6 > f5) { a1 = a6; a2 = a5; a3 = a4; f1 = f6; f2 = f5; f3 = f4; }
							else { a2 = a6; a3 = a5; f2 = f6; f3 = f5; }
						}
					}
				}
			}
			else {
				if (a4 - a2 > smallgapa * (a3 - a2)) {
					if (f2 >= f4) { a1 = a2; a2 = a4; f1 = f2; f2 = f4; }
					else { a3 = a4; f3 = f4; }
				}
				else {
					if (f4 > f2) {
						a5 = (a1 + a2) / 2; f5 = this.fun_ls(a5, x0, p);
						if (f5 > f2) { a1 = a5; a3 = a4; f1 = f5; f3 = f4; }
						else { a3 = a2; a2 = a5; f3 = f2; f2 = f5; }
					}
					else {
						a5 = (a3 + a4) / 2; f5 = this.fun_ls(a5, x0, p);
						if (f5 >= f4) { a1 = a2; a2 = a4; a3 = a5; f1 = f2; f2 = f4; f3 = f5; }
						else {
							a6 = (a3 + a5) / 2; f6 = this.fun_ls(a6, x0, p);
							if (f6 > f5) { a1 = a4; a2 = a5; a3 = a6; f1 = f4; f2 = f5; f3 = f6; }
							else { a1 = a5; a2 = a6; f1 = f5; f2 = f6; }
						}
					}
				}
			}
		}

		if (f2 > f0 && f4 > f0) a4 = 0;
		if (f2 <= f4) a4 = a2

		return a4;
	}

	H_end(x0, x1, f0, f1, e1, e2) {
		/*   Himmelblau termination rule.   return 1 for stop, 0 otherwise.
		*/
		let r = VectorUtils.norm(x0);
		if (r < e2) r = 1;
		r *= e1;
		if (VectorUtils.norm(VectorUtils.sub(x1, x0)) >= r)
			return 0;

		r = Math.abs(f0);
		if (r < e2) r = 1;
		r *= e1;
		if (Math.abs(f1 - f0) >= r)
			return 0;

		return 1;
	}

	/**
	* n-variate minimization with bounds using the BFGS algorithm
	* @param {number[]} x params
	* @param {[number, number][]} xb param bounds
	* @returns {number[]} minimized params
	*/
	ming2(x, xb) {
		let maxround = 100, maxstep = 8;
		let Ngoodtimes = 2, goodtimes = 0;
		let small = 1e-10, e = 1e-8, sizep0 = 0;
		let f0, ss;
		let n = x.length

		//xmark[n]: 0 for inside space; -1 for lower boundary; 1 for upper boundary.
		let xmark = x.map((_, i) => x[i] < xb[i][0] ? -1 : x[i] > xb[i][1] ? 1 : 0)
		x = x.map((v, i) => v < xb[i][0] ? xb[i][0] : v > xb[i][1] ? xb[i][1] : v)

		let xfree = xmark.map((x, i) => x == 0 ? -1 : i).filter(x => x != -1)
		f0 = this.lfun2dSdN(x);

		let x0 = [...x];
		let sizep = 999;

		let g0 = this.gradientB(x0, f0, xmark, sizep);

		let HH = MatrixUtils.diag(n);
		for (let i of xfree) HH[i][i] = 0;

		for (let Iround = 0; Iround < maxround; Iround++) {
			let p = HH.map(x => -VectorUtils.inner(x, g0))

			sizep0 = sizep;
			sizep = VectorUtils.norm(p)
			let am = maxstep, h
			for (let i = 0; i < n; i++)   /* max step length */
				if (p[i] > 0 && (xb[i][1] - x0[i]) / p[i] < am)
					am = (xb[i][1] - x0[i]) / p[i];
				else if (p[i] < 0 && (xb[i][0] - x0[i]) / p[i] < am)
					am = (xb[i][0] - x0[i]) / p[i];

			if (Iround == 0) {
				h = Math.abs(2 * f0 * 0.01 / VectorUtils.inner(g0, p));  /* check this?? */
				if (h > am / 2000) h = am / 2000
			}
			else {
				h = VectorUtils.norm(ss) / sizep;
				if (h < am / 500) h = am / 500
			}
			if (h < 1e-5) h = 1e-5;
			if (h > am / 5) h = am / 5;
			let alpha = this.LineSearch2(f0, x0, p, h, am, e);

			x = x0.map((_, i) => x0[i] + alpha * p[i])
			let f = this.lfun2dSdN(x)

			if (Iround == 0 || sizep < sizep0 || (sizep < .001 && sizep0 < .001))
				goodtimes++;
			else
				goodtimes = 0;
			if ((n == 1 || goodtimes >= Ngoodtimes) && sizep < 0.001 && f0 - f < 0.001 && this.H_end(x0, x, f0, f, e, e))
				break;

			let g = this.gradientB(x, f, xmark, sizep);
			/* modify the working set */
			for (let i = 0; i < n; i++) { /* add constraints, reduce H */

				if (xmark[i]) continue;

				if (Math.abs(x[i] - xb[i][0]) < 1e-6 && g[i] > 0)
					xmark[i] = -1;
				else if (Math.abs(x[i] - xb[i][1]) < 1e-6 && g[i] < 0)
					xmark[i] = 1;

				if (xmark[i] == 0) continue;

				for (let j in HH) {
					HH[i][j] = HH[j][i] = 0
				}
				xfree.push(i)

			}

			let [w, it] = g.reduce((a, _, i) => xmark[i] * g[i] > a[0] ? [xmark[i] * g[i], i] : a, [0, 0]) /* delete a constraint, enlarge H */
			if (w > 10 * sizep / (n - xfree.length)) {
				HH[it][it] = 1
				xfree = xfree.filter(x => x != it)

				xmark[it] = 0;
			}

			let yy = g.map((_, i) => xfree.includes(i) ? 0 : g[i] - g0[i])
			ss = x.map((_, i) => x[i] - x0[i])

			let zz = HH.map((_, i) => VectorUtils.inner(yy, HH[i]))
			let uu = VectorUtils.inner(yy, zz)

			f0 = f; g0 = [...g]; x0 = [...x]

			let vv = VectorUtils.inner(yy, ss)
			if (Math.abs(vv) < small) {
				HH = MatrixUtils.diag(n);
				for (let i of xfree) HH[i][i] = 0;
				continue;
			}

			for (let i = 0; i < n; i++) {
				for (let j = 0; j < n; j++) {
					HH[i][j] += ((1 + uu / vv) * ss[i] * ss[j] - zz[i] * ss[j] - ss[i] * zz[j]) / vv;
				}
			}
		}//end of for(Iround...)
		return x;
	}


	HouseholderRealSym(a, n) {

		let m, k, j, i;
		let scale, hh, h, g, f;
		let d = [], e = []

		for (i = n - 1; i >= 1; i--) {
			m = i - 1;
			h = scale = 0;
			if (m > 0) {
				for (k = 0; k <= m; k++)
					scale += Math.abs(a[i * n + k]);
				if (scale == 0)
					e[i] = a[i * n + m];
				else {
					for (k = 0; k <= m; k++) {
						a[i * n + k] /= scale;
						h += a[i * n + k] * a[i * n + k];
					}
					f = a[i * n + m];
					g = (f >= 0 ? -Math.sqrt(h) : Math.sqrt(h));
					e[i] = scale * g;
					h -= f * g;
					a[i * n + m] = f - g;
					f = 0;
					for (j = 0; j <= m; j++) {
						a[j * n + i] = a[i * n + j] / h;
						g = 0;
						for (k = 0; k <= j; k++)
							g += a[j * n + k] * a[i * n + k];
						for (k = j + 1; k <= m; k++)
							g += a[k * n + j] * a[i * n + k];
						e[j] = g / h;
						f += e[j] * a[i * n + j];
					}
					hh = f / (h * 2);
					for (j = 0; j <= m; j++) {
						f = a[i * n + j];
						e[j] = g = e[j] - hh * f;
						for (k = 0; k <= j; k++)
							a[j * n + k] -= (f * e[k] + g * a[i * n + k]);
					}
				}
			}
			else
				e[i] = a[i * n + m];
			d[i] = h;
		}
		d[0] = e[0] = 0;

		/* Get eigenvectors */
		for (i = 0; i < n; i++) {
			m = i - 1;
			if (d[i]) {
				for (j = 0; j <= m; j++) {
					g = 0;
					for (k = 0; k <= m; k++)
						g += a[i * n + k] * a[k * n + j];
					for (k = 0; k <= m; k++)
						a[k * n + j] -= g * a[k * n + i];
				}
			}
			d[i] = a[i * n + i];
			a[i * n + i] = 1;
			for (j = 0; j <= m; j++) a[j * n + i] = a[i * n + j] = 0;
		}

		return [d, e]
	}

	EigenTridagQLImplicit(d, e, n, z) {
		let m, j, iter, niter = 30, i, k;
		let s, r, p, g, f, dd, c, b, aa, bb;

		for (i = 1; i < n; i++) e[i - 1] = e[i]; e[n - 1] = 0;
		for (j = 0; j < n; j++) {
			iter = 0;
			do {
				for (m = j; m < n - 1; m++) {
					dd = Math.abs(d[m]) + Math.abs(d[m + 1]);
					if (Math.abs(e[m]) + dd == dd) break;  /* ??? */
				}
				if (m != j) {
					if (iter++ == niter) {
						break;
					}
					g = (d[j + 1] - d[j]) / (2 * e[j]);

					/* r=pythag(g,1); */

					if ((aa = Math.abs(g)) > 1) r = aa * Math.sqrt(1 + 1 / (g * g));
					else r = Math.sqrt(1 + g * g);

					g = d[m] - d[j] + e[j] / (g + (g >= 0 ? Math.abs(r) : -Math.abs(r)));
					s = c = 1;
					p = 0;
					for (i = m - 1; i >= j; i--) {
						f = s * e[i];
						b = c * e[i];

						/*  r=pythag(f,g);  */
						aa = Math.abs(f); bb = Math.abs(g);
						if (aa > bb) { bb /= aa; r = aa * Math.sqrt(1 + bb * bb); }
						else if (bb == 0) r = 0;
						else { aa /= bb; r = bb * Math.sqrt(1 + aa * aa); }

						e[i + 1] = r;
						if (r == 0) {
							d[i + 1] -= p;
							e[m] = 0;
							break;
						}
						s = f / r;
						c = g / r;
						g = d[i + 1] - p;
						r = (d[i] - g) * s + 2 * c * b;
						d[i + 1] = g + (p = s * r);
						g = c * r - b;
						for (k = 0; k < n; k++) {
							f = z[k * n + i + 1];
							z[k * n + i + 1] = s * z[k * n + i] + c * f;
							z[k * n + i] = c * z[k * n + i] - s * f;
						}
					}
					if (r == 0 && i >= j) continue;
					d[j] -= p; e[j] = g; e[m] = 0;
				}
			} while (m != j);
		}
		return [d, e, z];
	}
	EigenSort(d, U, n) {
		let k, j, i;
		let p;

		for (i = 0; i < n - 1; i++) {
			p = d[k = i];
			for (j = i + 1; j < n; j++)
				if (d[j] >= p) p = d[k = j];
			if (k != i) {
				d[k] = d[i];
				d[i] = p;
				for (j = 0; j < n; j++) {
					p = U[j * n + i];
					U[j * n + i] = U[j * n + k];
					U[j * n + k] = p;
				}
			}
		}
		return [d, U]
	}

	eigenRealSym(A, n) {
		let [Root, work] = this.HouseholderRealSym(A, n);
		[Root, work, A] = this.EigenTridagQLImplicit(Root, work, n, A);
		[Root, A] = this.EigenSort(Root, A, n);

		return [A, Root, work];
	}


	eigenQREV(Q, pi) {
		let i, j, Root = [], U = [], V = [], n = pi.length;

		let pi_sqrt = pi.map(x => Math.sqrt(x));

		/* store in U the symmetrical matrix S = sqrt(D) * Q * sqrt(-D) */
		for (i = 0; i < n; i++)
			for (j = 0, U[i * n + i] = Q[i * n + i]; j < i; j++)
				U[i * n + j] = U[j * n + i] = (Q[i * n + j] * pi_sqrt[i] / pi_sqrt[j]);

		[U, Root, V] = this.eigenRealSym(U, n);
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)  V[i * n + j] = U[j * n + i] * pi_sqrt[j];
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)  U[i * n + j] /= pi_sqrt[i];

		return [Root, U, V];
	}

	/**
	* Parse substitution rates according to the given model
	* @param {string} model
	* @param {string} kappa compressed substitution rate
	* @returns {number[]} `KAPPA` full substitution rate
	*/
	parseSubRates(model, kappa) {
		if (model == "JC" || model == "F81") {
			return [1, 1, 1, 1, 1, 1]
		}
		else if (model == "K2P" || model == "HKY") {//one para.
			return [kappa[0], kappa[0], 1, 1, 1, 1]
		}
		else if (model == "TNEF" || model == "TN") {//two para.
			return [kappa[0], kappa[1], 1, 1, 1, 1]
		}
		else if (model == "K3P" || model == "K3PUF") {//two para.
			return [kappa[0], kappa[0], kappa[1], kappa[1], 1, 1]
		}
		else if (model == "TIMEF" || model == "TIM") {//three para.
			return [kappa[0], kappa[1], kappa[2], kappa[2], 1, 1]
		}
		else if (model == "TVMEF" || model == "TVM") {//four para.
			return [kappa[0], kappa[0], kappa[1], kappa[2], kappa[3], 1]
		}
		return [kappa[0], kappa[1], kappa[2], kappa[3], kappa[4], 1]
	}

	EigenQc(omega, kappas, pi) {
		let n = pi.length;
		let pijQij;

		let ra0 = 0, rs = 0, ra = 0
		//Construct Q: transition probability matrix 64*64
		for (let i = 0; i < n; i++) {
			let codon1 = this.decodeCodon(i);
			for (let j = 0; j < i; j++) {
				let codon2 = this.decodeCodon(j);
				let diff = [0, 1, 2].filter(i => codon1[i] != codon2[i]);
				if (diff.length != 1) continue;
				let pos = diff[0]

				let ka = kappas[{
					"01": 0, "23": 1, "02": 2, "13": 3, "03": 4, "12": 5
				}[this.getSubstitutionType(codon1[pos], codon2[pos])]]

				//probability
				pijQij = 2 * pi[i] * ka * pi[j];

				if (this.isSynonyms(codon1, codon2)) {//synonymous
					rs += pijQij;
				}
				else {//nonsynonymous
					ra0 += pijQij;
					ra += pijQij * omega;
				}
			}
		}
		return [rs, ra, rs, ra0];
	}


	/**
	* Calculate decomposed transition probability matrix
	* @param {string} omega estimated Ka/Ks
	* @param {string} kappas estimated substitution rate
	* @param {string} pi64 codon frequency
	* @returns {[number[][], number[], number[]]} `[Root, U, V]`, where `[U] [Q] [U^-1] = [Root]` and `Q` is the transition probability matrix
	*/
	GetEigenQ(omega, kappas, pi64) {
		let pi = pi64.filter(x => x > 0)
		let n = pi.length
		let Q = new Array(n * n).fill(0);

		//Parse substitution rates according to the given model
		let rs = 0, ra = 0;
		//Construct Q: transition probability matrix 64*64
		for (let i = 0; i < n; i++) {
			let codon1 = this.decodeCodon(this.FROM61[i]);

			for (let j = 0; j < i; j++) {
				let codon2 = this.decodeCodon(this.FROM61[j]);
				let diff = [0, 1, 2].filter(i => codon1[i] != codon2[i]);
				//consider only one difference between two codons
				if (diff.length != 1) continue;
				let pos = diff[0]

				let n1 = codon1[pos],
					n2 = codon2[pos]
				let ka = kappas[{
					"01": 0, "23": 1, "02": 2, "13": 3, "03": 4, "12": 5
				}[this.getSubstitutionType(n1, n2)]]

				Q[i * n + j] = ka * pi[j];
				Q[j * n + i] = ka * pi[i];

				//probability
				let pijQij = 2 * pi[i] * Q[i * n + j];

				if (this.isSynonyms(codon1, codon2)) //synonymous
					rs += pijQij;
				else {//nonsynonymous
					Q[i * n + j] *= omega;
					Q[j * n + i] *= omega;
					ra += pijQij * omega;
				}

			} /* for (j) */
		}    /* for (i) */

		let mr = rs + ra;

		for (let i = 0; i < n; i++)
			Q[i * n + i] = -VectorUtils.sum(Q.slice(i * n, i * n + n));

		let [Root, U, V] = this.eigenQREV(Q, pi);
		Root = Root.map(x => x / mr);
		return [Root, U, V];
	}

	getSubstitutionType(nucl1, nucl2) {
		if (nucl1 != +nucl1) nucl1 = this.DNAencode[nucl1]
		if (nucl2 != +nucl2) nucl2 = this.DNAencode[nucl2]
		return [nucl1, nucl2].sort().join("")
	}

	/**
	* Likelihood function for calculating dS and dN
	* @param {number[]} param `[t, omega, ...kappas]` params
	* @returns {number} `lnL` maximum-likelihood score
	*/
	lfun2dSdN([t, omega, ...kappas]) {
		let n = this.FROM61.length;
		kappas = this.parseSubRates(this.model, kappas);

		let [Root, U, V] = this.GetEigenQ(omega, kappas, this.pi64);

		let expt = Root.map(r => Math.exp(t * r));

		let lnL1 = 0
		for (let c1 in this.fp) {
			let z0 = this.FROM64[this.encodeCodon(c1)]
			for (let c2 in this.fp[c1]) {
				let z1 = this.FROM64[this.encodeCodon(c2)]
				let fh = 0;
				for (let k = 0; k < n; k++) {
					fh += U[z0 * n + k] * expt[k] * V[k * n + z1];
				}
				fh *= this.pi64[this.encodeCodon(c1)];
				lnL1 -= Math.log(fh) * this.fp[c1][c2];
			}
		}
		return lnL1;
	}

	/**
	* calculate Ka/Ks
	* @param {string} seq1 encoded DNA sequence1
	* @param {string} seq2 encoded DNA sequence2
	* @returns {{Ka: number, Ks: number, S: number, N: number, Sd: number, Nd: number, t: number, method: string, KAPPA: number[], lnL: number, AICc: number}}
	*/
	calc(seq1, seq2) {
		let seq_length = seq1.length;
		let snp = VectorUtils.seq(seq_length).filter(i => seq1[i] != seq2[i]).length;

		let fp = {}
		let arr1 = seq1.match(/.{3}/g),
			arr2 = seq2.match(/.{3}/g)

		for (let h in arr1) {
			let codon1 = arr1[h], codon2 = arr2[h]
			if (codon1 < codon2) [codon1, codon2] = [codon2, codon1]
			if (!fp[codon1]) fp[codon1] = {}
			if (!fp[codon1][codon2]) fp[codon1][codon2] = 0
			fp[codon1][codon2]++
		}
		this.fp = fp

		if (["JC", "K2P", "TNEF", "K3P", "TIMEF", "TVMEF", "SYM"].includes(this.model)) {
			this.pi64 = new Array(64).fill(0)
			this.FROM61.forEach(i => this.pi64[i] = 1 / this.FROM61.length);
		} else {
			this.pi64 = this.getFreqency(seq1 + seq2)
		}

		/* initial values and bounds */
		let xbound = [
			[0.5 * snp / seq1.length, 3], //t > snp/length
			[0.001, 50],
			...new Array(this.nkappa).fill(0).map(_ => [0.01, 30])
		]

		//divergence time t, Ka/Ks, kappas
		let t = snp / seq1.length, omega = 0.25,
			kappas = new Array(this.nkappa).fill(1)

		let x = this.ming2([t, omega, ...kappas], xbound);

		[t, omega, ...kappas] = x;

		let lnL = -this.lfun2dSdN(x)

		//Parse substitution rates according to the given model
		let KAPPA = this.parseSubRates(this.model, kappas);

		let [rs, ra, rs0, ra0] = this.EigenQc(omega, KAPPA, this.pi64);
		let S = rs0 / (rs0 + ra0) * seq_length, N = seq_length - S,
			Sd = snp * rs / (rs + ra), Nd = snp - Sd,
			Ks = t / 3 * rs / (rs + ra) * (rs0 + ra0) / rs0,
			Ka = t / 3 * ra / (rs + ra) * (rs0 + ra0) / ra0

		//AICc = -2log(lnL) + 2K + 2K(K+1)/(n-K-1), K=parameters' number, n=sample size
		let AICc = -2 * lnL + 2 * (this.nkappa + 2) * (seq_length / 3) / ((seq_length / 3) - (this.nkappa + 2) - 1);

		return {
			Ka, Ks, S, N, Sd, Nd,
			t: (S * Ks + N * Ka) / (S + N),
			lnL, AICc, KAPPA,
			method: this.method
		};
	}

}
