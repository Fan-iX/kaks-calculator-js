/*********************************************************
* Filename: MYN.js
* Abstract: Definition of Modified YN00 (MYN) class.

* Version: js-2.0
* Author: Fanix
* Date: December.11, 2022

* Adapted from: KaKs_Calculator/MYN.cpp and KaKs_Calculator/MYN.h
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.

* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009

* Version: 1.0
* Author: Zhang Zhang (zhang.zhang@yale.edu)
* Date: Dec.30, 2005

* Modified Version: 2.0.1
* Modified Author: Kristian K Ullrich
* Modified Date: April.29, 2020

  References: 
	Zhang Zhang, Jun Li, Jun Yu. (2006) Computing Ka and Ks 
	with a consideration of unequal transitional substitutions. 
	BMC Evolutionary Biology, 6:44.
**********************************************************/

class MYN extends YN00 {
	constructor() {
		super();
		this.name = "MYN";
	}

	/* Get the two kappas between purines and between pyrimidines */
	GetKappa(seq1, seq2) {

		let i, j, k, h, pos, c = [], aa = [], b = [[], []], nondeg, fourdeg, by = [16, 4, 1];
		let kappatc_TN = [], kappaag_TN = [], kappa_TN = [];
		let S = [], wk = [], pi4 = [];
		let T1, T2, V;//proportions of transitional differences between purines and between
		let kdefault = 2, nullValue = 0.0;

		let F = [new Array(XSIZE).fill(0), new Array(XSIZE).fill(0)]

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


		for (k = 0; k < 2; k++) {  /* two kinds of sites */

			S[k] = this.sumArray(F[k], 16);
			if (S[k] <= 0) {
				wk[k] = 0;
				continue;
			}
			for (j = 0; j < 16; j++) {
				F[k][j] /= S[k];
			}

			//Transitions between purines
			T1 = 2 * F[k][2 * DNASIZE + 3];
			//Transitions between pyrimidines
			T2 = 2 * F[k][0 * DNASIZE + 1];
			//Tranversions
			V = 1 - T1 - T2 - (F[k][0 * DNASIZE + 0] + F[k][1 * DNASIZE + 1] + F[k][2 * DNASIZE + 2] + F[k][3 * DNASIZE + 3]);

			//pi[]: the sum probabilty of T, C, A, G, respectively
			for (j = 0; j < 4; j++) {
				pi4[j] = this.sumArray(F[k].slice(j * 4), 4);
			}

			[kappatc_TN[k], kappaag_TN[k]] = this.CorrectKappaTN93(S[k], T1, T2, V, pi4);
			wk[k] = ((kappatc_TN[k] > 0 && kappaag_TN[k] > 0) ? S[k] : 0);

			//R = (πTπCκ1 + πAπGκ2)/(πYπR), kappa = 2R in PAML's DOC
			kappa_TN[k] = 2 * (kappatc_TN[k] * pi4[0] * pi4[1] + kappaag_TN[k] * pi4[2] * pi4[3]) / ((pi4[0] + pi4[1]) * (pi4[2] + pi4[3]));

		}

		if (wk[0] + wk[1] == 0) {
			this.kappatc = this.kappaag = this.kappa = kdefault;
		}
		else {
			this.kappatc = (kappatc_TN[0] * wk[0] + kappatc_TN[1] * wk[1]) / (wk[0] + wk[1]);
			this.kappaag = (kappaag_TN[0] * wk[0] + kappaag_TN[1] * wk[1]) / (wk[0] + wk[1]);
			this.kappa = (kappa_TN[0] * wk[0] + kappa_TN[1] * wk[1]) / (wk[0] + wk[1]);
		}

		this.KAPPA[0] = this.kappatc;
		this.KAPPA[1] = this.kappaag;

		return 0;
	}


	/* Calculate transition probability matrix(64*64) */
	GetPMatCodon(kappa, omega) {

		let i, j, k, ndiff, pos = 0, from = new Array(3).fill(0), to = new Array(3).fill(0);
		let mr;
		let c = [];
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
					if ((from[pos] + to[pos] - 1) * (from[pos] + to[pos] - 5) == 0) {
						if (from[pos] + to[pos] == 1)
							PMatrix[i * CODON + j] *= this.kappatc;//T<->C
						else
							PMatrix[i * CODON + j] *= this.kappaag;//A<->G
					}

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
			//PMatrix[i*CODON+i] =- sumArray(PMatrix+i*CODON,CODON);
			PMatrix[i * CODON + i] -= this.sumArray(PMatrix.slice(i * CODON), CODON);
			//The sum of transition probability of main diagnoal elements
			mr -= this.pi[i] * PMatrix[i * CODON + i];
		}

		//calculate exp(PMatrix*t)
		[Root, U, V] = this.eigenQREV(PMatrix, this.pi, this.pi_sqrt, CODON, this.npi0)
		for (i = 0; i < CODON; i++)
			Root[i] /= mr;
		PMatrix = this.PMatUVRoot(this.t, CODON, U, V, Root);

		return PMatrix;
	}



	/* Correct kappas */
	CorrectKappaTN93(n, P1, P2, Q, pi4) {

		let failTN93;
		let tc, ag, Y, R, a1 = 0, a2 = 0, b = 0, A, B, C;
		let Qsmall = Math.min(1e-10, 0.1 / n), default_kappa = 2, maxkappa = 99;

		let kappatc_TN93 = -1, kappaag_TN93 = -1;
		failTN93 = 0;

		Y = pi4[0] + pi4[1];
		R = pi4[2] + pi4[3];

		tc = pi4[0] * pi4[1];
		ag = pi4[2] * pi4[3];


		if ((P1 + P2 + Q) > 1 || (P1 + P2) < -1e-10 || Q < -1e-10 || Math.abs(Y + R - 1) > 1e-8) {
			return [kappatc_TN93, kappaag_TN93];
		}

		if (Q < Qsmall)
			failTN93 = 1;
		else if (Y <= 0 || R <= 0 || (tc <= 0 && ag <= 0))
			failTN93 = 1;
		else {	//TN93 for multiple substitutions
			A = tc / Y + ag / R; B = tc + ag; C = Y * R;
			a1 = 1 - R * P1 / (2 * ag) - Q / (2 * R);
			a2 = 1 - Y * P2 / (2 * tc) - Q / (2 * Y);
			b = 1 - Q / (2 * C);
			if (a1 < 0 || a2 < 0 || b < 0) {
				failTN93 = 1;
			}
			else {
				if ((this.GAMMA == 20) || (this.GAMMA == -1)) {
					this.name = "GMYN";
				}

				if (this.GAMMA == 20) {
					a1 = Math.pow(a1, -1.0 / this.GAMMA);
					a2 = Math.pow(a2, -1.0 / this.GAMMA);
					b = Math.pow(b, -1.0 / this.GAMMA);
					//Kappa
					kappaag_TN93 = (R + Y * b - a1) / (R - R * b);
					kappatc_TN93 = (Y + R * b - a2) / (Y - Y * b);
				}
				else {
					a1 = Math.log(a1);
					a2 = Math.log(a2);
					b = Math.log(b);
					//Kappa
					kappaag_TN93 = (Y * b - a1) / (-R * b);
					kappatc_TN93 = (R * b - a2) / (-Y * b);
				}
			}
		}

		if (failTN93) {	//Fail to correct kappa
			kappatc_TN93 = kappaag_TN93 = default_kappa;
		}

		return [kappatc_TN93, kappaag_TN93]
	}

	/* Correct Ka and Ks */
	CorrectKaksTN93(n, P1, P2, Q, pi4, kaks, SEkaks) {

		let failTN93;
		let tc, ag, Y, R, a1, a2, b, A, B, C;
		let Qsmall = 1e-10;

		a1 = a2 = b = failTN93 = 0;

		Y = pi4[0] + pi4[1];
		R = pi4[2] + pi4[3];

		tc = pi4[0] * pi4[1];
		ag = pi4[2] * pi4[3];

		if (P1 + P2 + Q > 1 || Math.abs(Y + R - 1) > Qsmall || Y <= 0 || R <= 0 || (tc <= 0 && ag <= 0)) {
			failTN93 = 1;
		}
		else {	//TN93 for multiple substitutions
			A = tc / Y + ag / R; B = tc + ag; C = Y * R;
			a1 = 1 - R * P1 / (2 * ag) - Q / (2 * R);
			a2 = 1 - Y * P2 / (2 * tc) - Q / (2 * Y);
			b = 1 - Q / (2 * C);
			if (a1 < 0 || a2 < 0 || b < 0) {
				failTN93 = 1;
			}
			else {
				if ((this.GAMMA == 20) || (this.GAMMA == -1)) {
					this.name = "GMYN";
				}

				if (this.GAMMA == 20) {
					a1 = Math.pow(a1, -1.0 / this.GAMMA);
					a2 = Math.pow(a2, -1.0 / this.GAMMA);
					b = Math.pow(b, -1.0 / this.GAMMA);
					//Ka or Ks
					kaks = (ag * a1 / R) + tc * a2 / Y + (C - ag * Y / R - tc * R / Y) * b - ag - tc - Y * R;
					kaks *= 2 * this.GAMMA;
					let cc1 = a1 / (1 - R * P1 / (2 * ag) - Q / (2 * R));
					let cc2 = a2 / (1 - Y * P2 / (2 * tc) - Q / (2 * Y));
					let cc3 = ag / square(R) * cc1;
					cc3 += tc / square(Y) * cc2;
					cc3 += ((square(pi4[2]) + square(pi4[3])) / (2 * square(R)) + (square(pi4[0]) + square(pi4[1])) / (2 * square(Y))) * b / (1 - Q / (2 * C));
					SEkaks = (square(cc1) * P1 + square(cc2) * P2 + square(cc3) * Q - square(cc1 * P1 + cc2 * P2 + cc3 * Q)) / n;
				}

				else {
					a1 = Math.log(a1);
					a2 = Math.log(a2);
					b = Math.log(b);
					//Ka or Ks
					kaks = (-2 * ag * a1 / R) + (-2 * tc * a2 / Y) + (-2 * (C - ag * Y / R - tc * R / Y) * b);

					let cc1 = 2 * ag * R / (2 * ag * R - R * R * P1 - ag * Q);
					let cc2 = 2 * tc * Y / (2 * tc * Y - Y * Y * P2 - tc * Q);
					let cc3 = 2 * ag * ag / (R * (2 * ag * R - R * R * P1 - ag * Q));
					cc3 += 2 * tc * tc / (Y * (2 * tc * Y - Y * Y * P2 - tc * Q));
					cc3 += (R * R * (Y * Y - 2 * tc) + Y * Y * (R * R - 2 * ag)) / (2 * R * R * Y * Y - R * Y * Q);
					SEkaks = (square(cc1) * P1 + square(cc2) * P2 + square(cc3) * Q - square(cc1 * P1 + cc2 * P2 + cc3 * Q)) / n;
				}
			}
		}

		if (failTN93 == 1) {	//Use YN00's correction for Ka, Ks
			[Qsmall, kaks, SEkaks] = this.DistanceF84(n, P1 + P2, Q, pi4, Qsmall, kaks, SEkaks);
		}

		return [kaks, SEkaks];
	}



	/* Count differences, considering different transitional pathways between purines and between pyrimidines */
	CountDiffs(seq1, seq2, PMatrix) {
		let h, i1, i2, i, k, transi, c = [], ct = [0, 0], by = [16, 4, 1], aa = [];
		let dmark = [], step = [], b = [[], []], bt1 = [], bt2 = [];
		let ndiff, npath, nstop, sts1path = [], sts2path = [], stvpath = [], nts1path = [], nts2path = [], ntvpath = [];
		let sts1, sts2, stv, nts1, nts2, ntv; /* syn ts & tv, nonsyn ts & tv for 2 codons */
		let ppath = [], sump, p;

		let Sdts1 = 0, Sdts2 = 0, Sdtv = 0, Ndts1 = 0, Ndts2 = 0, Ndtv = 0;
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
			sts1 = sts2 = stv = nts1 = nts2 = ntv = 0;
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
				//transi=(transi==1 || transi==5);
				if (aa[0] == aa[1]) {
					if (transi == 5)
						sts1++;
					else if (transi == 1)
						sts2++;
					else
						stv++;
				}
				else {
					if (transi == 5)
						nts1++;
					else if (transi == 1)
						nts2++;
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
						step[0] = k / 2 |0;
						step[1] = k % 2;
						if (step[0] <= step[1])
							step[1]++;
						step[2] = 3 - step[0] - step[1];
					}//end of set the step[]

					for (i1 = 0; i1 < 3; i1++)
						bt1[i1] = bt2[i1] = b[0][i1];

					sts1path[k] = sts2path[k] = stvpath[k] = nts1path[k] = nts2path[k] = ntvpath[k] = 0;

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

						//ts & tr when syn & nonsyn in path k
						if (aa[0] == aa[1]) {
							if (transi == 5)
								sts1path[k]++;
							else if (transi == 1)
								sts2path[k]++;
							else
								stvpath[k]++;
						}
						else {
							if (transi == 5)
								nts1path[k]++;
							else if (transi == 1)
								nts2path[k]++;
							else
								ntvpath[k]++;
						}

						for (i2 = 0; i2 < 3; i2++)
							bt1[i2] = bt2[i2];
					}

				}  /* for(k,npath) */
				if (npath == nstop) {  /* all paths through stop codons */
					if (ndiff == 2) {
						nts1 = 0.25;
						nts2 = 0.25;
						ntv = 1.5;
					}
					else {
						nts1 = 0.25;
						nts2 = 0.25;
						ntv = 2.5;
					}
				}
				else {
					//sum probabilty of all path
					sump = this.sumArray(ppath, npath);
					if (sump > 1e-20) {
						for (k = 0; k < npath; k++) { //p: the probabilty of path k
							p = ppath[k] / sump;
							sts1 += sts1path[k] * p; sts2 += sts2path[k] * p; stv += stvpath[k] * p;
							nts1 += nts1path[k] * p; nts2 += nts2path[k] * p; ntv += ntvpath[k] * p;
						}
					}
				}

			}//end of if(ndiff)

			Sdts1 += sts1; Sdts2 += sts2; Sdtv += stv;
			Ndts1 += nts1; Ndts2 += nts2; Ndtv += ntv;
		}//end of for(h)

		return [Sdts1, Sdts2, Sdtv, Ndts1, Ndts2, Ndtv];
	}

	DistanceYN00(seq1, seq2) {

		let j, ir, nround = 100, status = 1;
		let fbS = new Array(4).fill(0), fbN = new Array(4).fill(0),
			fbSt, fbNt;
		let St, Nt, Sdts1, Sdts2, Sdtv, Ndts1, Ndts2, Ndtv;
		let w0 = 0, S0 = 0, N0 = 0, dS0 = 0, dN0 = 0, accu = 5e-8, minomega = 1e-5, maxomega = 99;
		let PMatrix = new Array(CODON * CODON).fill(0);
		let dS = 0.1, dN = 0.1, SEdS = NA, SEdN = NA;

		//initial values for t and omega(Ka/Ks)
		this.t = 0.09;
		this.omega = .5;
		this.S = this.N = 0.0;

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

		//Iterative loop
		for (ir = 0; ir < nround; ir++) {   /* iteration */

			//Get transition probability matrix from one codon to another
			PMatrix = this.GetPMatCodon(this.kappa, this.omega);

			//Count differences
			[Sdts1, Sdts2, Sdtv, Ndts1, Ndts2, Ndtv] = this.CountDiffs(seq1, seq2, PMatrix);

			//Synonymous(Sd) and nonsynonymous(Nd) differences
			this.Sd = Sdts1 + Sdts2 + Sdtv;
			this.Nd = Ndts1 + Ndts2 + Ndtv;

			//Seldom happen
			if (this.Sd > this.S) {
				Sdts1 *= (this.S / this.Sd);
				Sdts2 *= (this.S / this.Sd);
				Sdtv *= (this.S / this.Sd);
			}
			if (this.Nd > this.N) {
				Ndts1 *= (this.N / this.Nd);
				Ndts2 *= (this.N / this.Nd);
				Ndtv *= (this.N / this.Nd);
			}

			//Ks
			[dS, SEdS] = this.CorrectKaksTN93(this.S, Sdts1 / this.S, Sdts2 / this.S, Sdtv / this.S, fbS, dS, SEdS);
			//Ka
			[dN, SEdN] = this.CorrectKaksTN93(this.N, Ndts1 / this.N, Ndts2 / this.N, Ndtv / this.N, fbN, dN, SEdN);


			status = -1;

			if (dS < 1e-9) {
				status = -1;
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

		if (ir == nround)
			status = -2;

		return [dS, dN, SEdS, SEdN];
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
					if (k + b[j] == 1 || k + b[j] == 5) {//transition
						if (k + b[j] == 1) r *= this.kappatc;
						else r *= this.kappaag;	//(k+b[j]==5)
					}

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
}
