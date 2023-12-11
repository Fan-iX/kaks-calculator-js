/************************************************************
* Filename: GY94.js
* Abstract: Definition of GY94 class.

* Version: js-2.0
* Author: Fanix
* Date: December.11, 2022

* Adapted from: KaKs_Calculator/GY94.cpp and KaKs_Calculator/GY94.h
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.

* Modified Version: 2.0.1
* Modified Author: Kristian K Ullrich
* Modified Date: April.29, 2020

* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009

* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Oct.2, 2005

* Note: Source codes are taken from codeml.c in PAML.

  References:
  Goldman N, Yang Z  (1994)  A codon-based model of nucleotide
  substitution for protein-coding DNA sequences. Mol. Biol. 
  Evol. 11:725-736.
*************************************************************/

let GeneticCode =
	[[13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, -1, 17,
		10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
		9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
		19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7], /* 0:universal */

	[13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17,
		10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
		9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, -1, -1,
		19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7], /* 1:vertebrate mt.*/

	[13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17,
		16, 16, 16, 16, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
		9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
		19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7], /* 2:yeast mt. */

	[13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17,
		10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
		9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
		19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7], /* 3:mold mt. */

	[13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17,
		10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
		9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 15, 15,
		19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7], /* 4:invertebrate mt. */

	[13, 13, 10, 10, 15, 15, 15, 15, 18, 18, 5, 5, 4, 4, -1, 17,
		10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
		9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
		19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7], /* 5:ciliate nuclear*/

	[13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17,
		10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
		9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 2, 11, 15, 15, 15, 15,
		19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7], /* 6:echinoderm mt.*/

	[13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 4, 17,
		10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
		9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
		19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7], /* 7:euplotid mt. */

	[13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, -1, 17,
		10, 10, 10, 15, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
		9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
		19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7],
	/* 8:alternative yeast nu.*/

	[13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17,
		10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
		9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 7, 7,
		19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7], /* 9:ascidian mt. */

	[13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, 5, 4, 4, -1, 17,
		10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
		9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
		19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7], /* 10:blepharisma nu.*/

	[1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4,
		5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8,
		9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12,
		13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16] /* 11:Ziheng's regular code */
	];                                         /* GeneticCode[icode][#codon] */


/********************************************
* Function: construction function
* Input Parameter: a candidate model
* Output: 
* Return Value:

* Note: (Equal, Unequal: substitution rates)
JC, F81:    rTC==rAG =rTA==rCG==rTG==rCA
K2P, HKY:   rTC==rAG!=rTA==rCG==rTG==rCA
TNEF,TN:    rTC!=rAG!=rTA==rCG==rTG==rCA
K3P, K3PUF: rTC==rAG!=rTA==rCG!=rTG==rCA
TIMEF, TIM:	rTC!=rAG!=rTA==rCG!=rTG==rCA
TVMEF, TVM: rTC==rAG!=rTA!=rCG!=rTG!=rCA
SYM, GTR: 	rTC!=rAG!=rTA!=rCG!=rTG!=rCA
*********************************************/

//Constructor
class GY94 extends Base {
	constructor(NulModel) {
		super()
		this.FROM61 = new Array(CODON).fill(0)
		this.FROM64 = new Array(CODON).fill(0)
		this.name = "GY-" + NulModel;
		this.Iround = 0;
		this.SIZEp = 0;
		this.Small_Diff = 1e-6;
		this.w_rndu = 123456757;
		this.com = {
			z: [],
			ns: 2,
			ls: 0, ngene: 0, npatt: 0,
			icode: 0, ncode: 0,
			np: 0, nkappa: 0, sspace: 0,
			fpatt: new Array(CODON * CODON).fill(0),
			kappa: 0,
			omage: 0,
			pi: new Array(CODON).fill(0),
			KAPPA: new Array(8).fill(0)
		};
		this.lnL = 0.0;

		this.com.icode = genetic_code - 1;
		if (this.com.icode > 11) this.com.icode = 0;
		this.Nsensecodon = this.getNumNonsense(this.com.icode);
		this.com.ncode = 64 - this.Nsensecodon;

		//NulModel is set among a set of candidate models.
		this.model = NulModel;

		if (this.model == "JC" || this.model == "F81") this.com.nkappa = 0;
		else if (this.model == "K2P" || this.model == "HKY") this.com.nkappa = 1;
		else if (this.model == "TNEF" || this.model == "TN") this.com.nkappa = 2;
		else if (this.model == "K3P" || this.model == "K3PUF") this.com.nkappa = 2;
		else if (this.model == "TIMEF" || this.model == "TIM") this.com.nkappa = 3;
		else if (this.model == "TVMEF" || this.model == "TVM") this.com.nkappa = 4;
		else if (this.model == "SYM" || this.model == "GTR") this.com.nkappa = 5;

		//parameters' number: plus another two parameters (Ka/Ks and t)
		this.com.np = 2 + this.com.nkappa;
	}


	rndu() {
		this.w_rndu = (this.w_rndu * 69069 + 1) % 2 ** 32;
		return this.w_rndu * 2 ** -32;
	}


	/* x[i]=x0[i] + t*p[i] */
	fun_ls(t, x0, p, x, n) {
		for (let i = 0; i < n; i++) x[i] = x0[i] + t * p[i];
		return (this.lfun2dSdN(x, n));
	}



	PatternWeight() {
		let fpatti, lst = this.com.ls, h, ht, j, k = -1;
		let gap = 3;
		let nb = 1;
		let zt = [], b1 = [0, 0]; /* b[][] data at site h */
		let nc = 65;

		fpatti = new Array(lst).fill(0);

		for (j = 0; j < this.com.ns; j++) {
			zt[j] = new Array(lst).fill(0);
		}

		this.com.npatt = 0;

		for (h = 0; h < this.com.ls; h++) {

			b1[0] = this.com.z[0][h];
			b1[1] = this.com.z[1][h];

			//console.log(b1[1]);

			for (ht = 0; ht < this.com.npatt; ht++) {
				if (b1[0] == zt[0][ht] && b1[1] == zt[1][ht])
					break;
			}
			if (ht == this.com.npatt) {
				zt[0][this.com.npatt] = b1[0];
				zt[1][this.com.npatt] = b1[1];

				ht = this.com.npatt++;
			}

			fpatti[ht]++;
		}


		for (h = 0; h < this.com.npatt; h++) {
			this.com.fpatt[h] = fpatti[h];
		}

		for (j = 0; j < this.com.ns; j++) {
			this.com.z = zt
		}

		return 0;
	}

	preProcess(seq1, seq2) {
		let i;
		this.com.kappa = 2;
		this.com.omega = 0.4;

		this.setmark_61_64();

		this.com.ls = seq1.length;

		for (i = 0, this.snp = 0; i < this.com.ls; i++)
			if (seq1[i] != seq2[i]) this.snp++;

		this.com.z[0] = seq1;
		this.com.z[1] = seq2;

		this.com.ls /= 3;

		this.EncodeSeqs();
		this.PatternWeight();

		return 0;
	}

	setmark_61_64() {
		let i, code = GeneticCode[this.com.icode];
		//int c[3],aa0,aa1, by[3]={16,4,1};
		//double nSilent, nStop, nRepl;

		this.Nsensecodon = 0;
		for (i = 0; i < 64; i++) {
			if (code[i] == -1) {
				this.FROM64[i] = -1;
			}
			else {
				this.FROM61[this.Nsensecodon] = i;
				this.FROM64[i] = this.Nsensecodon++;
			}
		}
		this.com.ncode = this.Nsensecodon;
		return 0;
	}

	EncodeSeqs() {
		let j, h, k, b = new Array(3).fill(0);

		//T C A G to 0 1 2 3
		for (j = 0; j < this.com.ns; j++) {
			this.com.z[j] = this.com.z[j].split("").map(x => this.convertChar(x));
		}

		//encode to codon 0-64
		for (j = 0; j < this.com.ns; j++) {
			for (h = 0; h < this.com.ls; h++) {
				b[0] = this.com.z[j][h * 3]; b[1] = this.com.z[j][h * 3 + 1]; b[2] = this.com.z[j][h * 3 + 2];
				k = b[0] * 16 + b[1] * 4 + b[2];
				this.com.z[j][h] = this.FROM64[k];
			}
		}
	}

	GetCodonFreqs(pi) {
		let n = this.com.ncode, i, j, ic, b = new Array(3).fill(0);
		let fb3x4 = new Array(12).fill(0), fb4 = new Array(4).fill(0);
		let flag;

		for (i = 0, flag = new Array(CODON).fill(0); i < this.str1.length; i += 3) {
			j = this.getID(this.str1.substr(i, 3));
			if (this.getAminoAcid(j) != '!') flag[j] = 1;

			j = this.getID(this.str2.substr(i, 3));
			if (this.getAminoAcid(j) != '!') flag[j] = 1;
		}

		//Whether sequences are long enough
		if (this.sumArray(flag, CODON) == this.com.ncode) {
			return pi;
		}

		//Codon frequency is estimated from nucleotide frequency(fb3x4) at the three positions
		for (i = 0; i < n; i++) {
			ic = this.FROM61[i];
			b[0] = ic / 16 | 0;
			b[1] = (ic / 4 | 0) % 4;
			b[2] = ic % 4;
			for (j = 0; j < 3; j++) {
				fb3x4[j * 4 + b[j]] += pi[i];
				fb4[b[j]] += pi[i] / 3;
			}
		}

		//use nul frequencies to get codon frequency
		for (i = 0; i < n; i++) {
			ic = this.FROM61[i]; b[0] = ic / 16 | 0; b[1] = (ic / 4 | 0) % 4; b[2] = ic % 4;
			pi[i] = fb3x4[b[0]] * fb3x4[4 + b[1]] * fb3x4[8 + b[2]];
		}

		this.scaleArray(1. / this.sumArray(pi, n), pi, n);

		return pi;
	}

	gradientB(n, x, f0, g, space, xmark) {
		/* f0=fun(x) is always provided.
		   xmark=0: central; 1: upper; -1: down
		*/
		let i, j;
		let eh0 = this.Small_Diff, eh;  /* eh0=1e-6 || 1e-7 */

		for (i = 0; i < n; i++) {
			eh = eh0 * (Math.abs(x[i]) + 1);
			if (xmark[i] == 0 && this.SIZEp < 1) {    //central 
				for (j = 0; j < n; j++) space[j] = space[j + n] = x[j];
				eh = Math.pow(eh, .67); space[i] -= eh; space[i + n] += eh;
				g[i] = (this.lfun2dSdN(space.slice(n), n) - this.lfun2dSdN(space, n)) / (eh * 2.0);
			}
			else {//forward or backward
				for (j = 0; j < n; j++) space[j + n] = x[j];
				if (xmark[i]) eh *= -xmark[i];
				space[i + n] += eh;
				g[i] = (this.lfun2dSdN(space.slice(n), n) - f0) / eh;
			}
		}
		return (0);
	}


	LineSearch2(f, x0, p, step, limit, e, space, n) {
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
		let ii = 0, maxround = 10, nsymb = 0;
		let x = space, factor = 4, small = 1e-6, smallgapa = 0.2;
		let a0, a1, a2, a3, a4 = -1, a5, a6, f0, f1, f2, f3, f4 = -1, f5, f6;

		a0 = a1 = 0; f1 = f0 = f;
		a2 = a0 + step;
		f2 = this.fun_ls(a2, x0, p, x, n);
		if (f2 > f1) {
			for (; ;) {
				step /= factor;
				if (step < small) return [0, f, x0, p, space];
				a3 = a2; f3 = f2;
				a2 = a0 + step; f2 = this.fun_ls(a2, x0, p, x, n);
				if (f2 <= f1) break;
			}
		}
		else {
			for (; ;) {
				step *= factor;
				if (step > limit) step = limit;
				a3 = a0 + step; f3 = this.fun_ls(a3, x0, p, x, n);

				//obtain a bracket
				if (f3 >= f2) {	//a1<a2<a3 and f1>f2<f3
					break;
				}

				a1 = a2; f1 = f2; a2 = a3; f2 = f3;
				if (step >= limit) {
					f = f3; return [a3, f, x0, p, space];
				}
			}
		}

		/* iteration by quadratic interpolation, fig 2.2.9-10 (pp 71-71) */
		for (ii = 0; ii < maxround; ii++) {	//a4 is the minimum from the parabola over (a1,a2,a3)

			//2.2.4
			a4 = (a2 - a3) * f1 + (a3 - a1) * f2 + (a1 - a2) * f3;
			if (Math.abs(a4) > 1e-100) {
				//2.2.3
				a4 = ((a2 * a2 - a3 * a3) * f1 + (a3 * a3 - a1 * a1) * f2 + (a1 * a1 - a2 * a2) * f3) / (2 * a4);
			}
			if (a4 > a3 || a4 < a1) {//out of range: whether a1<a4<a3
				a4 = (a1 + a2) / 2;
			}
			f4 = this.fun_ls(a4, x0, p, x, n);

			if (Math.abs(f2 - f4) < e * (1 + Math.abs(f2))) {
				break;
			}

			if (a4 <= a2) {    /* fig 2.2.10 */
				if (a2 - a4 > smallgapa * (a2 - a1)) {
					if (f4 <= f2) { a3 = a2; a2 = a4; f3 = f2; f2 = f4; }
					else { a1 = a4; f1 = f4; }
				}
				else {
					if (f4 > f2) {
						a5 = (a2 + a3) / 2; f5 = this.fun_ls(a5, x0, p, x, n);
						if (f5 > f2) { a1 = a4; a3 = a5; f1 = f4; f3 = f5; }
						else { a1 = a2; a2 = a5; f1 = f2; f2 = f5; }
					}
					else {
						a5 = (a1 + a4) / 2; f5 = this.fun_ls(a5, x0, p, x, n);
						if (f5 >= f4) { a3 = a2; a2 = a4; a1 = a5; f3 = f2; f2 = f4; f1 = f5; }
						else {
							a6 = (a1 + a5) / 2; f6 = this.fun_ls(a6, x0, p, x, n);
							if (f6 > f5) { a1 = a6; a2 = a5; a3 = a4; f1 = f6; f2 = f5; f3 = f4; }
							else { a2 = a6; a3 = a5; f2 = f6; f3 = f5; }
						}
					}
				}
			}
			else {                     /* fig 2.2.9 */
				if (a4 - a2 > smallgapa * (a3 - a2)) {
					if (f2 >= f4) { a1 = a2; a2 = a4; f1 = f2; f2 = f4; }
					else { a3 = a4; f3 = f4; }
				}
				else {
					if (f4 > f2) {
						a5 = (a1 + a2) / 2; f5 = this.fun_ls(a5, x0, p, x, n);
						if (f5 > f2) { a1 = a5; a3 = a4; f1 = f5; f3 = f4; }
						else { a3 = a2; a2 = a5; f3 = f2; f2 = f5; }
					}
					else {
						a5 = (a3 + a4) / 2; f5 = this.fun_ls(a5, x0, p, x, n);
						if (f5 >= f4) { a1 = a2; a2 = a4; a3 = a5; f1 = f2; f2 = f4; f3 = f5; }
						else {
							a6 = (a3 + a5) / 2; f6 = this.fun_ls(a6, x0, p, x, n);
							if (f6 > f5) { a1 = a4; a2 = a5; a3 = a6; f1 = f4; f2 = f5; f3 = f6; }
							else { a1 = a5; a2 = a6; f1 = f5; f2 = f6; }
						}
					}
				}
			}
		}

		if (f2 > f0 && f4 > f0) a4 = 0;
		if (f2 <= f4) { f = f2; a4 = a2; }
		else f = f4;

		return [a4, f, x0, p, space];
	}

	distance(x, y, n) {
		let i, t = 0;
		for (i = 0; i < n; i++) t += square(x[i] - y[i]);
		return Math.sqrt(t);
	}

	H_end(x0, x1, f0, f1, e1, e2, n) {
		/*   Himmelblau termination rule.   return 1 for stop, 0 otherwise.
		*/
		let r;

		if ((r = this.norm(x0, n)) < e2) r = 1;
		r *= e1;
		if (this.distance(x1, x0, n) >= r)
			return 0;

		r = Math.abs(f0);
		if (r < e2) r = 1;
		r *= e1;
		if (Math.abs(f1 - f0) >= r)
			return 0;

		return 1;
	}

	ming2(f, x, xb, e, n) {
		/* n-variate minimization with bounds using the BFGS algorithm
			g0[n] g[n] p[n] x0[n] y[n] s[n] z[n] H[n*n] C[n*n] tv[2*n]
			xmark[n],ix[n]
			Size of space should be (check carefully?)
			#define spaceming2(n) ((n)*((n)*2+9+2)*sizeof(double))
			nfree: # free variables
			xmark[i]=0 for inside space; -1 for lower boundary; 1 for upper boundary.
			x[] has initial values at input and returns the estimates in return.
			ix[i] specifies the i-th free parameter
		*/
		let i, j, i1, i2, it, maxround = 2000, fail = 0, xmark, ix, nfree;
		let Ngoodtimes = 2, goodtimes = 0;
		let small = 1e-6, sizep0 = 0;
		let f0, g0, g, p, x0, y, s, z, H, C, tv;
		let w, v, alpha, am, h, maxstep = 8;

		if (n == 0) return [f, x, xb];

		g0 = new Array(n).fill(0); g = new Array(n).fill(0); x0 = new Array(n).fill(0);
		y = new Array(n).fill(0); s = new Array(n).fill(0); z = new Array(n).fill(0); H = new Array(n * n).fill(0);
		C = new Array(n * n).fill(0); tv = new Array(2 * n).fill(0);

		xmark = new Array(n).fill(0);
		ix = new Array(n).fill(0);

		for (i = 0; i < n; i++) {
			xmark[i] = 0;
			ix[i] = i;
		}

		for (i = 0, nfree = 0; i < n; i++) {
			if (x[i] <= xb[i][0]) {
				x[i] = xb[i][0];
				xmark[i] = -1;
				continue;
			}
			if (x[i] >= xb[i][1]) {
				x[i] = xb[i][1];
				xmark[i] = 1;
				continue;
			}

			ix[nfree++] = i;
		}

		f0 = f = this.lfun2dSdN(x, n);

		x0 = x.map(x => x);
		this.SIZEp = 999;

		this.gradientB(n, x0, f0, g0, tv, xmark);

		this.initIdentityMatrix(H, nfree);

		for (this.Iround = 0; this.Iround < maxround; this.Iround++) {
			for (i = 0, p = new Array(n).fill(0); i < nfree; i++) {
				for (j = 0; j < nfree; j++) p[ix[i]] -= H[i * nfree + j] * g0[ix[j]];
			}

			sizep0 = this.SIZEp;
			this.SIZEp = this.norm(p, n);      /* check this */

			for (i = 0, am = maxstep; i < n; i++) {  /* max step length */
				if (p[i] > 0 && (xb[i][1] - x0[i]) / p[i] < am)
					am = (xb[i][1] - x0[i]) / p[i];
				else if (p[i] < 0 && (xb[i][0] - x0[i]) / p[i] < am)
					am = (xb[i][0] - x0[i]) / p[i];
			}

			if (this.Iround == 0) {
				h = Math.abs(2 * f0 * .01 / this.innerp(g0, p, n));  /* check this?? */
				h = min2(h, am / 2000);
			}
			else {
				h = this.norm(s, nfree) / this.SIZEp;
				h = max2(h, am / 500);
			}
			h = max2(h, 1e-5); h = min2(h, am / 5);
			f = f0;
			[alpha, f, x0, p, tv] = this.LineSearch2(f, x0, p, h, am, min2(1e-3, e), tv, n); /* n or nfree? */

			fail = 0;
			for (i = 0; i < n; i++)  x[i] = x0[i] + alpha * p[i];

			w = min2(2, e * 1000);
			if (e < 1e-4 && e > 1e-6)
				w = 0.01;

			if (this.Iround == 0 || this.SIZEp < sizep0 || (this.SIZEp < .001 && sizep0 < .001))
				goodtimes++;
			else
				goodtimes = 0;
			if ((n == 1 || goodtimes >= Ngoodtimes) && this.SIZEp < (e > 1e-5 ? .05 : .001) && (f0 - f < 0.001) && this.H_end(x0, x, f0, f, e, e, n))
				break;

			this.gradientB(n, x, f, g, tv, xmark);

			/* modify the working set */
			for (i = 0; i < n; i++) {         /* add constraints, reduce H */

				if (xmark[i])
					continue;

				if (Math.abs(x[i] - xb[i][0]) < 1e-6 && -g[i] < 0)
					xmark[i] = -1;
				else if (Math.abs(x[i] - xb[i][1]) < 1e-6 && -g[i] > 0)
					xmark[i] = 1;

				if (xmark[i] == 0)
					continue;

				C = H;
				for (it = 0; it < nfree; it++)
					if (ix[it] == i)
						break;
				for (i1 = it; i1 < nfree - 1; i1++)
					ix[i1] = ix[i1 + 1];
				for (i1 = 0, nfree--; i1 < nfree; i1++)
					for (i2 = 0; i2 < nfree; i2++)
						H[i1 * nfree + i2] = C[(i1 + (i1 >= it)) * (nfree + 1) + i2 + (i2 >= it)];
			}

			for (i = 0, it = 0, w = 0; i < n; i++) {  /* delete a constraint, enlarge H */
				if (xmark[i] == -1 && -g[i] > w) {
					it = i; w = -g[i];
				}
				else if (xmark[i] == 1 && -g[i] < -w) {
					it = i; w = g[i];
				}
			}

			if (w > 10 * this.SIZEp / nfree) {
				C = H;

				for (i1 = 0; i1 < nfree; i1++) {
					for (i2 = 0; i2 < nfree; i2++) {
						H[i1 * (nfree + 1) + i2] = C[i1 * nfree + i2];
					}
				}

				for (i1 = 0; i1 < nfree + 1; i1++) {
					H[i1 * (nfree + 1) + nfree] = H[nfree * (nfree + 1) + i1] = 0;
				}

				H[(nfree + 1) * (nfree + 1) - 1] = 1;
				xmark[it] = 0;
				ix[nfree++] = it;
			}

			for (i = 0, f0 = f; i < nfree; i++) {
				y[i] = g[ix[i]] - g0[ix[i]];
				s[i] = x[ix[i]] - x0[ix[i]];
			}
			for (i = 0; i < n; i++) {
				g0[i] = g[i];
				x0[i] = x[i];
			}

			//reduce loop
			for (i = 0, w = v = 0.; i < nfree; i++) {
				for (j = 0, z[i] = 0.; j < nfree; j++) {
					z[i] += H[i * nfree + j] * y[j];
				}
				w += y[i] * z[i];
				v += y[i] * s[i];
			}
			if (Math.abs(v) < small) {
				this.initIdentityMatrix(H, nfree);
				fail = 1;
				continue;
			}
			for (i = 0; i < nfree; i++) {
				for (j = 0; j < nfree; j++) {
					H[i * nfree + j] += ((1 + w / v) * s[i] * s[j] - z[i] * s[j] - s[i] * z[j]) / v;
				}
			}

		}//end of for(Iround...)

		f = this.lfun2dSdN(x, n);

		return [f, x, xb];
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

					g = d[m] - d[j] + e[j] / (g + SIGN(r, g));
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


	eigenQREV(Q, pi, n, spacesqrtpi) {

		let i, j, inew, jnew, nnew, Root = [], U = [], V = [];
		let pi_sqrt = spacesqrtpi, small = 1e-6;

		for (j = 0, nnew = 0; j < n; j++)
			if (pi[j] > small)
				pi_sqrt[nnew++] = Math.sqrt(pi[j]);

		/* store in U the symmetrical matrix S = sqrt(D) * Q * sqrt(-D) */

		if (nnew == n) {
			for (i = 0; i < n; i++)
				for (j = 0, U[i * n + i] = Q[i * n + i]; j < i; j++)
					U[i * n + j] = U[j * n + i] = (Q[i * n + j] * pi_sqrt[i] / pi_sqrt[j]);

			[U, Root, V] = this.eigenRealSym(U, n);
			for (i = 0; i < n; i++) for (j = 0; j < n; j++)  V[i * n + j] = U[j * n + i] * pi_sqrt[j];
			for (i = 0; i < n; i++) for (j = 0; j < n; j++)  U[i * n + j] /= pi_sqrt[i];
		}
		else {
			for (i = 0, inew = 0; i < n; i++) {
				if (pi[i] > small) {
					for (j = 0, jnew = 0; j < i; j++)
						if (pi[j] > small) {
							U[inew * nnew + jnew] = U[jnew * nnew + inew]
								= Q[i * n + j] * pi_sqrt[inew] / pi_sqrt[jnew];
							jnew++;
						}
					U[inew * nnew + inew] = Q[i * n + i];
					inew++;
				}
			}

			[U, Root, V] = this.eigenRealSym(U, nnew);

			for (i = n - 1, inew = nnew - 1; i >= 0; i--)   /* construct Root */
				Root[i] = (pi[i] > small ? Root[inew--] : 0);
			for (i = n - 1, inew = nnew - 1; i >= 0; i--) {  /* construct V */
				if (pi[i] > small) {
					for (j = n - 1, jnew = nnew - 1; j >= 0; j--)
						if (pi[j] > small) {
							V[i * n + j] = U[jnew * nnew + inew] * pi_sqrt[jnew];
							jnew--;
						}
						else
							V[i * n + j] = +(i == j);
					inew--;
				}
				else
					for (j = 0; j < n; j++)  V[i * n + j] = +(i == j);
			}
			for (i = n - 1, inew = nnew - 1; i >= 0; i--) {  /* construct U */
				if (pi[i] > small) {
					for (j = n - 1, jnew = nnew - 1; j >= 0; j--)
						if (pi[j] > small) {
							U[i * n + j] = U[inew * nnew + jnew] / pi_sqrt[inew];
							jnew--;
						}
						else
							U[i * n + j] = +(i == j);
					inew--;
				}
				else
					for (j = 0; j < n; j++)
						U[i * n + j] = + (i == j);
			}
		}

		return [Root, U, V];
	}

	/********************************************
	* Function: parseSubRates
	* Input Parameter: string, array of double
	* Output: Parse substitution rates according to the given model
	* Return Value: int
	*********************************************/
	parseSubRates(model, kappa) {
		/*
		JC, F81:    rTC==rAG =rTA==rCG==rTG==rCA
		K2P, HKY:   rTC==rAG!=rTA==rCG==rTG==rCA
		TNEF,TN:    rTC!=rAG!=rTA==rCG==rTG==rCA
		K3P, K3PUF: rTC==rAG!=rTA==rCG!=rTG==rCA
		TIMEF, TIM:	rTC!=rAG!=rTA==rCG!=rTG==rCA
		TVMEF, TVM: rTC==rAG!=rTA!=rCG!=rTG!=rCA
		SYM, GTR: 	rTC!=rAG!=rTA!=rCG!=rTG!=rCA
		*/

		kappa[5] = 1.0;		//Substitution rate between C and A

		if (this.model == "JC" || this.model == "F81") {
			//Q[i*n+j] = kappa[0] = 1;
			kappa[0] = kappa[1] = kappa[2] = kappa[3] = kappa[4] = kappa[5];
		}
		else if (this.model == "K2P" || this.model == "HKY") {//one para.
			//if ((b1+b2)==1 || (b1+b2)==5) Q[i*n+j] = kappa[0];
			kappa[1] = kappa[0];
			kappa[2] = kappa[3] = kappa[4] = kappa[5];
		}
		else if (this.model == "TNEF" || this.model == "TN") {//two para.
			//if ((b1+b2)==1)  Q[i*n+j] = kappa[0];
			//else if ((b1+b2)==5)  Q[i*n+j] = kappa[1];
			kappa[2] = kappa[3] = kappa[4] = kappa[5];
		}
		else if (this.model == "K3P" || this.model == "K3PUF") {//two para.
			//if ((b1+b2)==1 || (b1+b2)==5)  Q[i*n+j] = kappa[0];
			//else if ((b1+b2)==2 || (b1+b2)==4) Q[i*n+j] = kappa[1];
			kappa[4] = kappa[5];
			kappa[2] = kappa[3] = kappa[1];
			kappa[1] = kappa[0];
		}
		else if (this.model == "TIMEF" || this.model == "TIM") {//three para.
			//if ((b1+b2)==1) Q[i*n+j] = kappa[0];
			//else if ((b1+b2)==5)  Q[i*n+j] = kappa[1];
			//else if ((b1+b2)==2 || (b1+b2)==4) Q[i*n+j] = kappa[2];
			kappa[4] = kappa[5];
			kappa[3] = kappa[2];
		}
		else if (this.model == "TVMEF" || this.model == "TVM") {//four para.
			//if ((b1+b2)==1 || (b1+b2)==5)  Q[i*n+j] = kappa[0];
			//else if ((b1+b2)==2) Q[i*n+j] = kappa[1];
			//else if ((b1+b2)==4) Q[i*n+j] = kappa[2];
			//else if ((b1==0) && (b2==3)) Q[i*n+j] = kappa[3];

			//TVMEF, TVM: rTC==rAG!=rTA!=rCG!=rTG!=rCA
			//SYM, GTR:   rTC!=rAG!=rTA!=rCG!=rTG!=rCA
			kappa[4] = kappa[3];
			kappa[3] = kappa[2];
			kappa[2] = kappa[1];
			kappa[1] = kappa[0];
		}
		else if (this.model == "SYM" || this.model == "GTR") {//five para.
		}
		else {
			return 0;
		}

		return 1;
	}

	/* Construct Q: transition probability matrix 64*64 */
	EigenQc(getstats, blength, kappa, omega) {
		/* This contructs the rate matrix Q for codon substitution and get the eigen
		   values and vectors if getstats==0, or get statistics (dS & dN) if 
		   getstats==1.  
		   The routine is also called by Qcodon2aa for mechanistic amino acid 
		   substitution models.
		   Input parameters are kappa, omega and com.pi (or com.fb61).
		
		   Statistics calculated include S, dS & dN.
		   c0[0,1,2] and c[0,1,2] are rates for the 3 codon positions before and after 
		   selection.  c4 is for 4-fold rates.  ts[3] and tv[3] are transition/
		   transversion rates, not calculated.
		
		   Under NSsites or other site-class models, this function does not scale 
		   Q but calculates the Qfactor_NS.
		   DetailOutput() uses this function to calculate dS & dN; for each omega 
		   component, dS and dN are not scaled here but are scaled in DetailOutput()
		   after Qfactor_NS is calculated.
		
		   aaDist=FIT1 & FIT2:  ap,p*,av,v*, (and w0 for FIT2)
		   The argument omega is used only if the model assumes one omega.  For 
		   AAClasses, com.pomega is used instead.
		*/
		let n = this.Nsensecodon, i, j, k, ic1, ic2, aa1, aa2;
		let ndiff, pos = 0, from = [0, 0, 0], to = [0, 0, 0];
		let mr, rs0, ra0, rs, ra; /* rho's */
		let d4 = 0, d0 = [0, 0, 0], d = [0, 0, 0], ts = [0, 0, 0], tv = [0, 0, 0];  /* rates at positions and 4-fold sites */
		let pi = this.com.pi, w = -1, pijQij; let S, dS, dN, Root, U, V;
		let space_pisqrt = new Array(CODON).fill(0);

		let Q = new Array(n * n).fill(0);

		//Equal codon frequency
		if (this.model == "JC" || this.model == "K2P" || this.model == "TNEF" || this.model == "K3P" ||
			this.model == "TIMEF" || this.model == "TVMEF" || this.model == "SYM") {
			pi = this.com.pi = new Array(64).fill(1.0 / (64 - this.getNumNonsense(genetic_code)));
		}

		//Parse substitution rates according to the given model
		this.parseSubRates(this.model, this.com.KAPPA);

		//Construct Q: transition probability matrix 64*64
		for (i = 0, rs0 = ra0 = rs = ra = 0; i < n; i++) {

			//codon i
			ic1 = this.FROM61[i]; from[0] = ic1 / 16 | 0; from[1] = (ic1 / 4 | 0) % 4; from[2] = ic1 % 4;

			//codon j
			for (j = 0; j < i; j++) {
				ic2 = this.FROM61[j]; to[0] = ic2 / 16 | 0; to[1] = (ic2 / 4 | 0) % 4; to[2] = ic2 % 4;
				for (k = 0, ndiff = 0; k < 3; k++) {
					if (from[k] != to[k]) {
						ndiff++; pos = k;
					}
				}

				//consider only one difference between two codons
				if (ndiff != 1) continue;

				let b1 = min2(from[pos], to[pos]);
				let b2 = max2(from[pos], to[pos]);

				/*           01   23   02   13   03   12
				JC, F81:    rTC==rAG =rTA==rCG==rTG==rCA
				K2P, HKY:   rTC==rAG!=rTA==rCG==rTG==rCA
				TNEF,TN:    rTC!=rAG!=rTA==rCG==rTG==rCA
				K3P, K3PUF: rTC==rAG!=rTA==rCG!=rTG==rCA
				TIMEF, TIM:	rTC!=rAG!=rTA==rCG!=rTG==rCA
				TVMEF, TVM: rTC==rAG!=rTA!=rCG!=rTG!=rCA
				SYM, GTR: 	rTC!=rAG!=rTA!=rCG!=rTG!=rCA
				*/
				if (b1 == 0 && b2 == 1) Q[i * n + j] = kappa[0]; /* TC */
				else if (b1 == 2 && b2 == 3) Q[i * n + j] = kappa[1]; /* AG */
				else if (b1 == 0 && b2 == 2) Q[i * n + j] = kappa[2]; /* TA */
				else if (b1 == 1 && b2 == 3) Q[i * n + j] = kappa[3]; /* CG */
				else if (b1 == 0 && b2 == 3) Q[i * n + j] = kappa[4]; /* TG */
				else if (b1 == 1 && b2 == 2) Q[i * n + j] = kappa[5]; /* CA=1 */

				Q[j * n + i] = Q[i * n + j];
				Q[i * n + j] *= this.com.pi[j];
				Q[j * n + i] *= this.com.pi[i];

				//probability
				pijQij = 2 * pi[i] * Q[i * n + j];

				aa1 = GeneticCode[this.com.icode][ic1];
				aa2 = GeneticCode[this.com.icode][ic2];
				if (aa1 == aa2) {//synonymous
					rs += pijQij;
				}
				else {//nonsynonymous
					ra0 += pijQij;
					w = omega;

					Q[i * n + j] *= w; Q[j * n + i] *= w;
					ra += pijQij * w;
				}

			} /* for (j) */
		}    /* for (i) */

		mr = rs + ra;
		if (getstats) {
			rs0 = rs;
			w = (rs0 + ra0); rs0 /= w; ra0 /= w; S = rs0 * 3 * this.com.ls;
			if (blength >= 0) {  /* calculates dS & dN */
				if (blength == 0) dS = dN = 0;
				rs /= mr;
				ra /= mr;
				dS = blength * rs / (3 * rs0);
				dN = blength * ra / (3 * ra0);
			}
			return [S, dS, dN, Q];
		}
		else {
			for (i = 0; i < n; i++) Q[i * n + i] = -this.sumArray(Q.slice(i * n), n);

			[Root, U, V] = this.eigenQREV(Q, this.com.pi, n, space_pisqrt);

			for (i = 0; i < n; i++) Root[i] /= mr;
			return [Root, U, V, Q];
		}
	}

	/* Return maximum-likelihood score */
	lfun2dSdN(x, np) {
		/* likelihood function for calculating dS and dN between 2 sequences,
		   com.z[0] & com.z[1]:
				 prob(i,j) = PI_i * p(i,j,t)
		   
		   Data are clean and coded.
		   Transition probability pijt is calculated for observed patterns only.
		*/
		let n = this.com.ncode, h, k, ik, z0, z1;
		let fh, expt = [], lnL1 = 0;
		let pkappa = this.com.KAPPA;

		//cout<<name.c_str()<<": "<<com.ncode<<"\t"<<Nsensecodon<<endl;

		k = 1, ik = 0;
		this.com.kappa = x[k];
		for (ik = 0; ik < this.com.nkappa; ik++)
			pkappa[ik] = x[k++];

		this.com.omega = x[1 + this.com.nkappa];

		let [Root, U, V, PMat] = this.EigenQc(0, -1, pkappa, this.com.omega);

		//t = x[0],  exp(Qt)
		for (k = 0; k < n; k++) {
			expt[k] = Math.exp(x[0] * Root[k]);
		}

		//com.npatt = number of patterns
		for (h = 0; h < this.com.npatt; h++) {

			if (this.com.fpatt[h] < this.Small_Diff) {
				continue;
			}

			z0 = this.com.z[0][h];
			z1 = this.com.z[1][h];

			for (k = 0, fh = 0; k < n; k++) {
				fh += U[z0 * n + k] * expt[k] * V[k * n + z1];
			}

			fh *= this.com.pi[z0];

			lnL1 -= Math.log(fh) * this.com.fpatt[h];
		}

		return lnL1;
	}

	/* Main function to calculate Ka and Ks, called by "Run" */
	PairwiseCodon() {
		let pz0 = []
		let npatt0 = this.com.npatt;
		let fpatt0, ls0 = this.com.ls;
		let fp = new Array(CODON * CODON).fill(0);
		let n = this.com.ncode, j, k, h;
		let x = [.3, 1, .5, .5, .5, .5, .3, 0, 0, 0],
			xb = [[1e-6, 3], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]];
		let kappab = [.01, 30], omegab = [0.001, 50];
		let e = 1e-6, dS, dN, PMat;
		let pkappa = this.com.KAPPA;

		fpatt0 = new Array(npatt0 * 3);
		pz0 = this.com.z.map(x => x.map(y => y))

		fpatt0 = this.com.fpatt.map(x => x);

		//t > snp/length
		if ((this.snp / length) > 1e-6) xb[0][0] = 3 * this.snp / length;

		//com.nkappa = 2;
		for (j = 0; j < this.com.nkappa; j++) {
			xb[1 + j][0] = kappab[0];
			xb[1 + j][1] = kappab[1];
		}

		k = 1 + this.com.nkappa;
		xb[k][0] = omegab[0];
		xb[k][1] = omegab[1];

		//fp: codon -> codon   frequency
		fp = new Array(CODON * CODON).fill(0);
		for (h = 0; h < npatt0; h++) {
			j = max2(pz0[0][h], pz0[1][h]);
			k = min2(pz0[0][h], pz0[1][h]);
			fp[j * n + k] += fpatt0[h];
		}

		//fp -> com.fpatt
		for (j = 0, this.com.npatt = 0; j < n; j++) {
			for (k = 0; k < j + 1; k++)
				if (fp[j * n + k]) {
					this.com.z[0][this.com.npatt] = j;
					this.com.z[1][this.com.npatt] = k;
					this.com.fpatt[this.com.npatt++] = fp[j * n + k];
				}
		}

		//com.pi: sense codon's frequencies
		for (j = 0, this.com.pi = new Array(n).fill(0); j < this.com.npatt; j++) {
			this.com.pi[this.com.z[0][j]] += this.com.fpatt[j] / (2 * this.com.ls);
			this.com.pi[this.com.z[1][j]] += this.com.fpatt[j] / (2 * this.com.ls);
		}

		this.com.pi = this.GetCodonFreqs(this.com.pi);

		/* initial values and bounds */
		//divergence time t
		x[0] = -1;
		if (x[0] > 3)
			x[0] = 1.5 + this.rndu();
		if (x[0] < 1e-6)
			x[0] = .5 * this.rndu();

		//kappas
		for (j = 0; j < this.com.nkappa; j++)
			x[1 + j] = .2 + .4 * this.rndu();

		//Ka/Ks
		k = 1 + this.com.nkappa;
		x[k] = (3 * x[k] + 0.6 * this.rndu()) / 4;
		x[k] = max2(x[k], 0.01);
		x[k] = min2(x[k], 2);

		if ((this.snp / length) > 1e-6) x[0] = 3.0 * this.snp / length;

		[this.lnL, x, xb] = this.ming2(this.lnL, x, xb, e, this.com.np);

		[this.S, dS, dN, PMat] = this.EigenQc(1, x[0], pkappa, this.com.omega);

		this.Ka = dN;
		this.Ks = dS;
		this.N = this.com.ls * 3 - this.S;
		this.Sd = this.Ks * this.S;
		this.Nd = this.Ka * this.N;
		let scale = (this.Sd + this.Nd) / this.snp;
		this.Sd /= scale;
		this.Nd /= scale;

		this.lnL = -this.lnL;
		//AICc = -2log(lnL) + 2K + 2K(K+1)/(n-K-1), K=parameters' number, n=sample size
		this.AICc = -2 * this.lnL + 2. * (this.com.nkappa + 2) * (length / 3) / ((length / 3) - (this.com.nkappa + 2) - 1.);

		this.t = x[0] / 3;
		//kappa = this.com.kappa;

		this.parseSubRates(this.model, pkappa);
		this.KAPPA = pkappa

		k = this.com.np - 1;

		this.com.ls = ls0;
		for (k = 0; k < this.com.ns; k++) this.com.z[k] = pz0[k];
		this.com.npatt = npatt0;

		for (h = 0; h < npatt0; h++) this.com.fpatt[h] = fpatt0[h];

		return 0;
	}

	/* Main fuction for GY method */
	Run(seq1, seq2) {
		this.str1 = seq1;
		this.str2 = seq2;

		this.preProcess(seq1, seq2);

		this.PairwiseCodon();

		return this.parseOutput();
	}

}
