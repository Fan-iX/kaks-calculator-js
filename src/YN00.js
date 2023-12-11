/************************************************************
* Filename: YN00.js
* Abstract: Definition of YN00 class.

* Version: js-2.0
* Author: Fanix
* Date: December.11, 2022

* Adapted from: KaKs_Calculator/YN00.cpp and KaKs_Calculator/YN00.h
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.

* Modified Version: 2.0.1
* Modified Author: Kristian K Ullrich
* Modified Date: April.29, 2020
 
* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009

* Version: 1.0
* Author: Zhang Zhang (zhang.zhang@yale.edu)
* Date: Jan.31, 2005

  Note: Source codes are adapted from yn00.c in PAML.

  Reference:
  Yang Z, Nielsen R  (2000)  Estimating Synonymous and 
  Nonsynonymous Substitution Rates Under Realistic 
  Evolutionary Models. Mol Biol Evol 17:32-43.
*************************************************************/
CODONFREQ = 12
class YN00 extends Base {
	constructor() {
		super();
		this.name = "YN";
		this.f12pos = new Array(12).fill(0);
		this.pi = new Array(64).fill(0);
		this.pi_sqrt = new Array(64).fill(0);
		this.iteration = 1;
	}

	getFreqency(seq1, seq2) {
		let i, fstop = 0;

		//Get A,C,G,T frequency at three positions
		for (i = 0; i < seq1.length; i++) {
			this.f12pos[(i % 3) * 4 + this.convertChar(seq1[i])]++;
			this.f12pos[(i % 3) * 4 + this.convertChar(seq2[i])]++;
		}
		for (i = 0; i < CODONFREQ; i++)
			this.f12pos[i] /= (seq1.length + seq2.length) / 3;

		//Get 64 amino acid probability
		for (i = 0; i < CODON; i++) {
			this.pi[i] = this.f12pos[i / 16 | 0] * this.f12pos[4 + (i % 16) / 4 | 0] * this.f12pos[8 + i % 4];
			if (this.getAminoAcid(i) == '!') {
				fstop += this.pi[i];
				this.pi[i] = 0;
			}
		}

		//Scale the sum of pi[] to 1
		for (i = 0; i < CODON; i++)
			this.pi[i] /= (1.0 - fstop);

		if (Math.abs(1 - this.sumArray(this.pi, CODON)) > 1e-6)
			console.warn("Warning: error in get codon freqency.");

		for (i = 0, this.npi0 = 0; i < CODON; i++)
			if (this.pi[i])
				this.pi_sqrt[this.npi0++] = Math.sqrt(this.pi[i]);

		this.npi0 = CODON - this.npi0;
	}

	/* Estimate kappa using the fourfold degenerate sites at third codon positions and nondegenerate sites */
	GetKappa(seq1, seq2) {

		let i, j, k, h, pos, c = [], aa = [], b = [[], []], nondeg, fourdeg, by = [16, 4, 1];
		let ka = [0, 0], S = [], wk = [], T, V, pi4 = [];
		let kdefault = 2, nullValue = 0.0, t;

		let F = [new Array(XSIZE).fill(0), new Array(XSIZE).fill(0)];

		//Get Pi[] of A,C,G,T
		for (h = 0; h < seq1.length; h += 3) {

			//c[]: amino acid(0--63)
			c[0] = this.getID(seq1.substr(h, 3));
			c[1] = this.getID(seq2.substr(h, 3));
			//aa[ ]: amino acid
			aa[0] = this.getAminoAcid(c[0]);
			aa[1] = this.getAminoAcid(c[1]);
			//b[][]: 0--3
			for (j = 0; j < 3; j++) {
				b[0][j] = this.convertChar(seq1[h + j]);
				b[1][j] = this.convertChar(seq2[h + j]);
			}

			//Find non-degenerate sites
			for (pos = 0; pos < 3; pos++) {
				for (k = 0, nondeg = 0; k < 2; k++) {
					for (i = 0; i < 4; i++) {
						if (i != b[k][pos])
							if (this.getAminoAcid(c[k] + (i - b[k][pos]) * by[pos]) == aa[k])
								break;
					}
					if (i == 4)
						nondeg++;
				}
				//F[0][]: 0-fold
				if (nondeg == 2) {
					F[0][b[0][pos] * 4 + b[1][pos]] += .5;
					F[0][b[1][pos] * 4 + b[0][pos]] += .5;
				}
			}

			//Find 4-fold degenerate sites at 3rd position
			for (k = 0, fourdeg = 0; k < 2; k++) {
				for (j = 0, i = c[k] - b[k][2]; j < 4; j++)
					if (j != b[k][2] && this.getAminoAcid(i + j) != aa[k])
						break;
				if (aa[0] == aa[1] && j == 4)
					fourdeg++;
			}
			//F[1][]: 4-fold
			if (fourdeg == 2) {
				F[1][b[0][2] * 4 + b[1][2]] += .5;
				F[1][b[1][2] * 4 + b[0][2]] += .5;
			}

		}//end of for(h)	

		for (k = 0; k < 2; k++) { /* two kinds of sites */

			S[k] = this.sumArray(F[k], 16);
			if (S[k] <= 0) {
				wk[k] = 0;
				continue;
			}
			for (j = 0; j < 16; j++)
				F[k][j] /= S[k];

			//Transition
			T = (F[k][0 * 4 + 1] + F[k][2 * 4 + 3]) * 2;
			//Tranversion
			V = 1 - T - (F[k][0 * 4 + 0] + F[k][1 * 4 + 1] + F[k][2 * 4 + 2] + F[k][3 * 4 + 3]);

			//pi[]: the sum probabilty of T, C, A, G, respectively
			for (j = 0; j < 4; j++) {
				pi4[j] = this.sumArray(F[k].slice(j * 4), 4);
			}

			//Correct kappa
			[ka[k], t, nullValue] = this.DistanceF84(S[k], T, V, pi4, ka[k], t, nullValue);
			wk[k] = (ka[k] > 0 ? S[k] : 0);
		}

		if (wk[0] + wk[1] == 0) {
			this.kappa = kdefault;
		}
		else {
			this.kappa = (ka[0] * wk[0] + ka[1] * wk[1]) / (wk[0] + wk[1]);
		}

		this.KAPPA[0] = this.KAPPA[1] = this.kappa;

		return 0;
	}


	/* Correct for multiple substitutions */
	DistanceF84(n, P, Q, pi4, k_HKY, t, SEt) {

		let failF84 = 0, failK80 = 0, failJC69 = 0;
		let tc, ag, Y, R, a = 0, b = 0, A, B, C, k_F84;
		let Qsmall = Math.min(1e-10, 0.1 / n), maxkappa = 2, maxt = 99;

		k_HKY = -1;

		Y = pi4[0] + pi4[1];
		R = pi4[2] + pi4[3];

		tc = pi4[0] * pi4[1];
		ag = pi4[2] * pi4[3];

		//wangdp added
		if (this.GAMMA == 4 || this.GAMMA == -1) {
			this.name = "GYN";
		}

		if (P + Q > 1) {
			t = maxt;
			k_HKY = 1;
			return [k_HKY, t, SEt];
		}
		if (P < -1e-10 || Q < -1e-10 || Math.abs(Y + R - 1) > 1e-8) {
			return [k_HKY, t, SEt];
		}

		//HKY85
		if (Q < Qsmall)
			failF84 = failK80 = 1;
		else if (Y <= 0 || R <= 0 || (tc <= 0 && ag <= 0))
			failF84 = 1;
		else {
			A = tc / Y + ag / R; B = tc + ag; C = Y * R;
			a = (2 * B + 2 * (tc * R / Y + ag * Y / R) * (1 - Q / (2 * C)) - P) / (2 * A);
			b = 1 - Q / (2 * C);
			if (a <= 0 || b <= 0)
				failF84 = 1;
		}

		if (!failF84) {

			if ((this.GAMMA == 4) || (this.GAMMA == 20)) //WANGDP
			{
				a = .5 * this.GAMMA * (Math.pow(a, -1.0 / this.GAMMA) - 1);
				b = .5 * this.GAMMA * (Math.pow(b, -1.0 / this.GAMMA) - 1);
			}
			else {
				a = -.5 * Math.log(a); b = -.5 * Math.log(b);
			}

			if (b <= 0)
				failF84 = 1;
			else {
				k_F84 = a / b - 1;
				t = 4 * b * (tc * (1 + k_F84 / Y) + ag * (1 + k_F84 / R) + C);
				k_HKY = (B + (tc / Y + ag / R) * k_F84) / B; /* k_F84=>k_HKY85 */

				//Standard errors
				a = A * C / (A * C - C * P / 2 - (A - B) * Q / 2);
				b = A * (A - B) / (A * C - C * P / 2 - (A - B) * Q / 2) - (A - B - C) / (C - Q / 2);
				SEt = Math.sqrt((a * a * P + b * b * Q - (a * P + b * Q) ** 2) / n);
			}
		}

		//K80
		if (failF84 && !failK80) {  /* try K80 */
			a = 1 - 2 * P - Q; b = 1 - 2 * Q;
			if (a <= 0 || b <= 0)
				failK80 = 1;
			else {
				a = -Math.log(a); b = -Math.log(b);
				if (b <= 0)
					failK80 = 1;
				else {
					k_HKY = (.5 * a - .25 * b) / (.25 * b);
					t = .5 * a + .25 * b;
				}
				if (SEt) {
					a = 1 / (1 - 2 * P - Q); b = (a + 1 / (1 - 2 * Q)) / 2;
					SEt = Math.sqrt((a * a * P + b * b * Q - (a * P + b * Q) ** 2) / n);
				}
			}
		}

		if (failK80) {/* try JC69 */
			if ((P += Q) >= .75) {
				failJC69 = 1;
				P = .75 * (n - 1.) / n;
			}
			t = -.75 * Math.log(1 - P * 4 / 3.);

			if (t > 99)
				t = maxt;
			if (SEt) {
				SEt = Math.sqrt(9 * P * (1 - P) / n) / (3 - 4 * P);
			}
		}

		if (k_HKY > 99)
			k_HKY = maxkappa;

		return [k_HKY, t, SEt];
	}

	DistanceYN00(seq1, seq2) {
		let j, k, ir, nround = 10;
		let fbS = new Array(4).fill(0), fbN = new Array(4).fill(0),
			fbSt, fbNt;
		let St, Nt, Sdts, Sdtv, Ndts, Ndtv, k_HKY = 0;
		let w0 = 0, dS0 = 0, dN0 = 0, accu = 5e-4, minomega = 1e-5, maxomega = 99;
		let PMatrix = new Array(CODON * CODON).fill(0);
		let dS = 0.1, dN = 0.1, SEKs = NA, SEKa = NA

		if (this.t == 0)
			this.t = .5;
		if (this.omega <= 0)
			this.omega = 1;
		for (k = 0; k < 4; k++)
			fbS[k] = fbN[k] = 0;

		//Count sites of sequence 1
		[St, Nt, fbSt, fbNt] = this.CountSites(seq1);
		this.S += St / 2;
		this.N += Nt / 2;
		for (j = 0; j < 4; j++) {
			fbS[j] += fbSt[j] / 2;
			fbN[j] += fbNt[j] / 2;
		}
		//Count sites of sequence 2
		[St, Nt, fbSt, fbNt] = this.CountSites(seq2);
		this.S += St / 2;
		this.N += Nt / 2;
		for (j = 0; j < 4; j++) {
			fbS[j] += fbSt[j] / 2;
			fbN[j] += fbNt[j] / 2;
		}

		if (this.t < 0.001 || this.t > 5)
			t = 0.5;
		if (this.omega < 0.01 || this.omega > 5)
			this.omega = .5;


		for (ir = 0; ir < (this.iteration ? nround : 1); ir++) {   /* iteration */
			if (this.iteration)
				PMatrix = this.GetPMatCodon(this.kappa, this.omega);
			else
				for (j = 0; j < CODON * CODON; j++) {
					PMatrix[j] = 1;
				}

			[Sdts, Sdtv, Ndts, Ndtv] = this.CountDiffs(seq1, seq2, PMatrix);

			this.Sd = Sdts + Sdtv;
			this.Nd = Ndts + Ndtv;

			//synonymous
			[k_HKY, dS, SEKs] = this.DistanceF84(this.S, Sdts / this.S, Sdtv / this.S, fbS, k_HKY, dS, SEKs);
			//nonsynonymous
			[k_HKY, dN, SEKa] = this.DistanceF84(this.N, Ndts / this.N, Ndtv / this.N, fbN, k_HKY, dN, SEKa);

			if (dS < 1e-9) {
				this.omega = maxomega;
			}
			else {
				this.omega = Math.max(minomega, dN / dS);
			}

			this.t = dS * 3 * this.S / (this.S + this.N) + dN * 3 * this.N / (this.S + this.N);

			if (Math.abs(dS - dS0) < accu && Math.abs(dN - dN0) < accu && Math.abs(this.omega - w0) < accu)
				break;

			dS0 = dS;
			dN0 = dN;
			w0 = this.omega;

		} //end of for(ir) */

		return [dS, dN, SEKs, SEKa];
	}

	//Count differences between two compared codons
	CountDiffs(seq1, seq2, PMatrix) {
		let h, i1, i2, i, k, transi, c = [], ct = [], by = [16, 4, 1], aa = [];
		let dmark = [], step = [], b = [[], []], bt1 = [], bt2 = [];
		let ndiff, npath, nstop, stspath = [], stvpath = [], ntspath = [], ntvpath = [];
		let sts, stv, nts, ntv; /* syn ts & tv, nonsyn ts & tv for 2 codons */
		let ppath = [], sump, p;
		let Sdts = 0, Sdtv = 0, Ndts = 0, Ndtv = 0

		this.snp = 0;
		for (h = 0; h < seq1.length; h += 3) {

			c[0] = this.getID(seq1.substr(h, 3))
			c[1] = this.getID(seq2.substr(h, 3));
			//Difference?
			if (c[0] == c[1])
				continue;

			for (i = 0; i < 2; i++) {
				b[i][0] = c[i] / 16 | 0;
				b[i][1] = (c[i] % 16) / 4 | 0;
				b[i][2] = c[i] % 4;
				aa[i] = this.getAminoAcid(c[i]);
			}

			//ndiff: differences of two codons
			ndiff = 0;
			sts = stv = nts = ntv = 0;
			//dmark[]: position of different codon 
			for (k = 0; k < 3; k++) {
				dmark[k] = -1;
				if (b[0][k] != b[1][k])
					dmark[ndiff++] = k;
			}

			this.snp += ndiff;

			npath = 1;
			if (ndiff > 1)
				npath = (ndiff == 2) ? 2 : 6;

			if (ndiff == 1) {
				transi = b[0][dmark[0]] + b[1][dmark[0]];
				transi = (transi == 1 || transi == 5);
				if (aa[0] == aa[1]) {
					if (transi)
						sts++;
					else
						stv++;
				}
				else {
					if (transi)
						nts++;
					else
						ntv++;
				}
			}
			else {   /* ndiff=2 or 3 */

				nstop = 0;
				for (k = 0; k < npath; k++) {

					//set the step[]
					for (i1 = 0; i1 < 3; i1++)
						step[i1] = -1;
					if (ndiff == 2) {
						step[0] = dmark[k];
						step[1] = dmark[1 - k];
					}
					else {
						step[0] = k / 2 | 0;
						step[1] = k % 2;
						if (step[0] <= step[1])
							step[1]++;
						step[2] = 3 - step[0] - step[1];
					}//end of set the step[]

					for (i1 = 0; i1 < 3; i1++)
						bt1[i1] = bt2[i1] = b[0][i1];

					stspath[k] = stvpath[k] = ntspath[k] = ntvpath[k] = 0;

					//ppath[]: probabilty of each path
					for (i1 = 0, ppath[k] = 1; i1 < ndiff; i1++) {
						bt2[step[i1]] = b[1][step[i1]];

						//ct[]: mutated codon's ID(0--63)
						for (i2 = 0, ct[0] = ct[1] = 0; i2 < 3; i2++) {
							ct[0] += bt1[i2] * by[i2];
							ct[1] += bt2[i2] * by[i2];
						}
						//ppath[k]: probabilty of path k
						ppath[k] *= PMatrix[ct[0] * CODON + ct[1]];
						for (i2 = 0; i2 < 2; i2++)
							aa[i2] = this.getAminoAcid(ct[i2]);

						if (aa[1] == '!') {
							nstop++;
							ppath[k] = 0;
							break;
						}

						transi = b[0][step[i1]] + b[1][step[i1]];
						transi = (transi == 1 || transi == 5);  /* transition? */

						//ts & tr when syn & nonsyn in path k
						if (aa[0] == aa[1]) {
							if (transi)
								stspath[k]++;
							else
								stvpath[k]++;
						}
						else {
							if (transi)
								ntspath[k]++;
							else
								ntvpath[k]++;
						}

						for (i2 = 0; i2 < 3; i2++)
							bt1[i2] = bt2[i2];
					}

				}  /* for(k,npath) */
				if (npath == nstop) {  /* all paths through stop codons */
					if (ndiff == 2) {
						nts = .5;
						ntv = 1.5;
					}
					else {
						nts = .5;
						ntv = 2.5;
					}
				}
				else {
					//sum probabilty of all path
					sump = this.sumArray(ppath, npath);
					if (sump > 1e-20) {
						for (k = 0; k < npath; k++) { //p: the probabilty of path k
							p = ppath[k] / sump;
							sts += stspath[k] * p; stv += stvpath[k] * p;
							nts += ntspath[k] * p; ntv += ntvpath[k] * p;
						}
					}
				}
			}//end of if(ndiff)
			Sdts += sts;
			Sdtv += stv;
			Ndts += nts;
			Ndtv += ntv;

		}//end of for(h)

		return [Sdts, Sdtv, Ndts, Ndtv];
	}

	/* Calculate transition probability matrix(64*64) */
	GetPMatCodon(kappa, omega) {
		/* Get PMat=exp(Q*t) for weighting pathways
			*/
		let i, j, k, ndiff, pos = 0, from = new Array(3).fill(0), to = new Array(3).fill(0);
		let mr;
		let c = [0, 0];
		let U = new Array(CODON * CODON).fill(0),
			V = new Array(CODON * CODON).fill(0),
			Root = new Array(CODON * CODON).fill(0),
			PMatrix = new Array(CODON * CODON).fill(0);
		for (i = 0; i < CODON; i++) {
			for (j = 0; j < i; j++) {

				//codon 'from'
				from[0] = i / 16 | 0; from[1] = (i / 4 | 0) % 4; from[2] = i % 4;
				//codon 'to'
				to[0] = j / 16 | 0; to[1] = (j / 4 | 0) % 4; to[2] = j % 4;
				//amino acid of 'from' and 'to'
				c[0] = this.getAminoAcid(i);
				c[1] = this.getAminoAcid(j);
				//stop codon
				if (c[0] == '!' || c[1] == '!')
					continue;

				//whether two codons only have one difference
				for (k = 0, ndiff = 0; k < 3; k++) {
					if (from[k] != to[k]) {
						ndiff++;
						pos = k;
					}
				}
				if (ndiff == 1) {
					//only have one difference
					PMatrix[i * CODON + j] = 1;
					//transition
					if ((from[pos] + to[pos] - 1) * (from[pos] + to[pos] - 5) == 0)
						PMatrix[i * CODON + j] *= kappa;

					//nonsynonymous
					if (c[0] != c[1])
						PMatrix[i * CODON + j] *= omega;

					//diagonal element is equal
					PMatrix[j * CODON + i] = PMatrix[i * CODON + j];
				}
			}
		}

		//PMatrix[](*Q): transition probability matrix
		for (i = 0; i < CODON; i++)
			for (j = 0; j < CODON; j++)
				PMatrix[i * CODON + j] *= this.pi[j];

		//scale the sum of PMat[][j](j=0-63) to zero
		for (i = 0, mr = 0; i < CODON; i++) {
			PMatrix[i * CODON + i] = -this.sumArray(PMatrix.slice(i * CODON), CODON);
			//The sum of transition probability of main diagnoal elements
			mr -= this.pi[i] * PMatrix[i * CODON + i];
		}

		//calculate exp(PMatrix*t)
		[Root, U, V] = this.eigenQREV(PMatrix, this.pi, this.pi_sqrt, CODON, this.npi0);
		for (i = 0; i < CODON; i++)
			Root[i] /= mr;
		PMatrix = this.PMatUVRoot(this.t, CODON, U, V, Root);

		return PMatrix;
	}



	/* P(t) = U * exp{Root*t} * V */
	PMatUVRoot(t, n, U, V, Root) {
		let i, j, k;
		let expt, uexpt;
		let smallp = 0;
		let P = new Array(n * n).fill(0)

		for (k = 0; k < n; k++)
			for (i = 0, expt = Math.exp(t * Root[k]); i < n; i++)
				for (j = 0, uexpt = U[i * n + k] * expt; j < n; j++)
					P[i * n + j] += uexpt * V[k * n + j];

		for (i = 0; i < n * n; i++)
			if (P[i] < smallp)
				P[i] = 0;

		return P;
	}


	/* Count the synonymous and nonsynonymous sites of two sequences */
	CountSites(seq) {
		let h, i, j, k, c = [], aa = [], b = [], by = [16, 4, 1];
		let r, S, N;

		let Stot = 0, Ntot = 0, fbS = [0, 0, 0, 0], fbN = [0, 0, 0, 0]

		for (h = 0; h < seq.length; h += 3) {

			//Get codon id and amino acid
			c[0] = this.getID(seq.substr(h, 3));
			aa[0] = this.getAminoAcid(c[0]);
			for (i = 0; i < 3; i++) {
				b[i] = this.convertChar(seq[h + i]);
			}

			for (j = 0, S = N = 0; j < 3; j++) {
				for (k = 0; k < 4; k++) {    /* b[j] changes to k */
					if (k == b[j])
						continue;
					//c[0] change at position j
					c[1] = c[0] + (k - b[j]) * by[j];
					aa[1] = this.getAminoAcid(c[1]);

					if (aa[1] == '!')
						continue;

					r = this.pi[c[1]];
					if (k + b[j] == 1 || k + b[j] == 5)	//transition
						r *= this.kappa;

					if (aa[0] == aa[1]) { //synonymous
						S += r;
						fbS[b[j]] += r; //syn probability of A,C,G,T					
					}
					else { //nonsynonymous
						N += r;
						fbN[b[j]] += r; //nonsyn probability of A,C,G,T					
					}
				}
			}
			Stot += S;
			Ntot += N;
		}

		//Scale Stot+Ntot to seq.length
		r = seq.length / (Stot + Ntot);
		Stot *= r;
		Ntot *= r;

		//get probablity of syn of four nul.
		r = this.sumArray(fbS, 4);
		for (k = 0; k < 4; k++)
			fbS[k] /= r;

		//get probablity of nonsyn of four nul.
		r = this.sumArray(fbN, 4);
		for (k = 0; k < 4; k++)
			fbN[k] /= r;

		return [Stot, Ntot, fbS, fbN];
	}


	Run(seq1, seq2) {
		this.t = 0.4;
		this.kappa = NA;
		this.omega = 1;
		this.Ks = this.Ka = 0.1;

		this.getFreqency(seq1, seq2);
		this.GetKappa(seq1, seq2);
		[this.Ks, this.Ka, this.SEKs, this.SEKa] = this.DistanceYN00(seq1, seq2);

		this.t = (this.S * this.Ks + this.N * this.Ka) / (this.S + this.N);

		return this.parseOutput();
	}



	//The following functions are used to calculate PMatrix by Taylor.

	eigenQREV(Q, pi, pi_sqrt, n, npi0) {
		/* 
		This finds the eigen solution of the rate matrix Q for a time-reversible 
		Markov process, using the algorithm for a real symmetric matrix.
		Rate matrix Q = S * diag{pi} = U * diag{Root} * V, 
		where S is symmetrical, all elements of pi are positive, and U*V = I.
		pi_sqrt[n-npi0] has to be calculated before calling this routine.
		
		  [U 0] [Q_0 0] [U^-1 0]    [Root  0]
		  [0 I] [0   0] [0    I]  = [0     0]
		  
			Ziheng Yang, 25 December 2001 (ref is CME/eigenQ.pdf)
			*/
		let i, j, inew, jnew, nnew = n - npi0, Root = [], U = [], V = [];

		//npi0 is the number of stop codons in selected genetic table

		if (this.npi0 == 0) {	//seldom occur

			//Set U[64*64]
			for (i = 0; i < n; i++) {
				for (j = 0, U[i * n + i] = Q[i * n + i]; j < i; j++)
					U[i * n + j] = U[j * n + i] = (Q[i * n + j] * pi_sqrt[i] / pi_sqrt[j]);
			}

			//Set U[64*64]
			[U, Root, V] = this.eigenRealSym(U, n);

			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
					V[i * n + j] = U[j * n + i] * pi_sqrt[j];

			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
					U[i * n + j] /= pi_sqrt[i];
		}
		else {
			for (i = 0, inew = 0; i < n; i++) {
				if (pi[i]) {
					for (j = 0, jnew = 0; j < i; j++)
						if (pi[j]) {
							U[inew * nnew + jnew] = U[jnew * nnew + inew] = Q[i * n + j] * pi_sqrt[inew] / pi_sqrt[jnew];
							jnew++;
						}
					U[inew * nnew + inew] = Q[i * n + i];
					inew++;
				}
			}
			[U, Root, V] = this.eigenRealSym(U, nnew);


			for (i = n - 1, inew = nnew - 1; i >= 0; i--)   /* construct Root */
				Root[i] = (pi[i] ? Root[inew--] : 0);

			for (i = n - 1, inew = nnew - 1; i >= 0; i--) {  /* construct V */
				if (pi[i]) {
					for (j = n - 1, jnew = nnew - 1; j >= 0; j--)
						if (pi[j]) {
							V[i * n + j] = U[jnew * nnew + inew] * pi_sqrt[jnew];
							jnew--;
						}
						else
							V[i * n + j] = +(i == j);
					inew--;
				}
				else
					for (j = 0; j < n; j++)
						V[i * n + j] = +(i == j);
			}

			for (i = n - 1, inew = nnew - 1; i >= 0; i--) {  /* construct U */
				if (pi[i]) {
					for (j = n - 1, jnew = nnew - 1; j >= 0; j--)
						if (pi[j]) {
							U[i * n + j] = U[inew * nnew + jnew] / pi_sqrt[inew];
							jnew--;
						}
						else
							U[i * n + j] = +(i == j);
					inew--;
				}
				else
					for (j = 0; j < n; j++)
						U[i * n + j] = +(i == j);
			}
		}

		return [Root, U, V];
	}

	HouseholderRealSym(a, n) {
		/* This uses HouseholderRealSym transformation to reduce a real symmetrical matrix 
		a[n*n] into a tridiagonal matrix represented by d and e.
		d[] is the diagonal (eigends), and e[] the off-diagonal.
			*/
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
		/* This finds the eigen solution of a tridiagonal matrix represented by d and e.  
		d[] is the diagonal (eigenvalues), e[] is the off-diagonal
		z[n*n]: as input should have the identity matrix to get the eigen solution of the 
		tridiagonal matrix, or the output from HouseholderRealSym() to get the 
		eigen solution to the original real symmetric matrix.
		z[n*n]: has the orthogonal matrix as output
		
		  Adapted from routine tqli in Numerical Recipes in C, with reference to
		  LAPACK fortran code.
		  Ziheng Yang, May 2001
			*/
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
		/* this sorts the eigen values d[] and rearrange the (right) eigen vectors U[]
			*/
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
		/* This finds the eigen solution of a real symmetrical matrix A[n*n].  In return, 
		A has the right vectors and Root has the eigenvalues. work[n] is the working space.
		The matrix is first reduced to a tridiagonal matrix using HouseholderRealSym(), 
		and then using the QL algorithm with implicit shifts.  
		
		  Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK
		  Ziheng Yang, 23 May 2001
			*/
		let [Root, work] = this.HouseholderRealSym(A, n);
		[Root, work, A] = this.EigenTridagQLImplicit(Root, work, n, A);
		[Root, A] = this.EigenSort(Root, A, n);

		return [A, Root, work];
	}
}
