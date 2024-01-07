/***********************************************************
* Filename: YN.js
* Abstract: Definition of YN and MYN class.

* Version: js-3.1
* Author: Fanix
* Date: December.24, 2023

* Re-implementation KaKs_Calculator/YN00 and KaKs_Calculator/MYN

  Reference:
  Yang Z, Nielsen R  (2000)  Estimating Synonymous and 
  Nonsynonymous Substitution Rates Under Realistic 
  Evolutionary Models. Mol Biol Evol 17:32-43.

  Yang Z  (1994)  Maximum likelihood phylogenetic estimation 
  from DNA sequences with variable rates over sites: 
  approximate methods. J. Mol. Evol. 39:306–314.
*************************************************************/
class YN extends Base {
	constructor(genetic_code = 1) {
		super(genetic_code);
		this.method = "YN";
	}

	/**
	* Test if two nucleotide is transition
	* @param {string} nucl1 encoded nucleotide
	* @param {string} nucl2 encoded nucleotide
	* @returns {boolean}
	*/
	isTransition(nucl1, nucl2) {
		return nucl1 == this.transite(nucl2)
	}

	/**
	* Count codon frequence
	* @param {string} seq encoded DNA sequence
	* @returns {number[]} `pi` codon frequence
	*/
	getFreqency(seq) {
		let pi = new Array(64).fill(0);

		//Get A,C,G,T frequency at three positions
		let f12pos = new Array(12).fill(0);
		(seq).split("").forEach((x, i) => f12pos[(i % 3) * 4 + +x]++)
		f12pos = f12pos.map(x => x / (seq.length / 3))

		// Get 64 amino acid probability
		for (let i of [0, 1, 2, 3])
			for (let j of [0, 1, 2, 3])
				for (let k of [0, 1, 2, 3])
					if (this.translation[i * 16 + j * 4 + k] != '!')
						pi[i * 16 + j * 4 + k] = f12pos[i] * f12pos[j + 4] * f12pos[k + 8];

		//Scale the sum of pi[] to 1
		let scale = VectorUtils.sum(pi)
		pi = pi.map(x => x / scale)
		return pi
	}

	/**
	* Estimate substitution rate (kappa) using the fourfold degenerate sites at third codon positions and nondegenerate sites
	* @param {string} seq1 encoded DNA sequence
	* @param {string} seq2 encoded DNA sequence
	* @returns {number} `kappa` substitution rate
	*/
	GetKappa(seq1, seq2) {
		let kdefault = 2;

		let F0 = new Array(4).fill(0).map(_ => new Array(4).fill(0)),
			F4 = new Array(4).fill(0).map(_ => new Array(4).fill(0));

		let arr1 = seq1.match(/.{3}/g),
			arr2 = seq2.match(/.{3}/g)
		//Get Pi[] of A,C,G,T

		for (let h in arr1) {
			let codon1 = arr1[h],
				codon2 = arr2[h];

			//Find non-degenerate sites
			for (let pos of [0, 1, 2])
				if (
					this.getMutatedCodonList(codon1, pos)
						.every(c => !this.isSynonyms(c, codon1)) &&
					this.getMutatedCodonList(codon2, pos)
						.every(c => !this.isSynonyms(c, codon2))
				) {
					F0[codon1[pos]][codon2[pos]] += .5;
					F0[codon2[pos]][codon1[pos]] += .5;
				}

			// Find 4-fold degenerate sites at 3rd position
			if (
				this.isSynonyms(codon1, codon2) &&
				this.getMutatedCodonList(codon1, 2)
					.every(c => this.isSynonyms(c, codon1)) &&
				this.getMutatedCodonList(codon2, 2)
					.every(c => this.isSynonyms(c, codon2))
			) {
				F4[codon1[2]][codon2[2]] += .5;
				F4[codon2[2]][codon1[2]] += .5;
			}
		}//end of for(h)	

		let freq2Kappa = (F) => {
			let S = VectorUtils.sum(F.flat());
			if (S < 0) return 0;
			F = F.map(x => x.map(y => y / S))
			//Transition
			let T = (F[0][1] + F[2][3]) * 2;
			//Tranversion
			let V = 1 - T - [0, 1, 2, 3].reduce((a, i) => a + F[i][i], 0);
			//pi[]: the sum probabilty of T, C, A, G, respectively
			let pi4 = [0, 1, 2, 3].map(i => VectorUtils.sum(F[i]));
			//Correct kappa
			let ka = this.CorrectKappaF84(T, V, pi4) || this.CorrectKappaK80(T, V);
			return [ka, ka > 0 ? S : 0];
		}

		let [kaF0, wkF0] = freq2Kappa(F0),
			[kaF4, wkF4] = freq2Kappa(F4);

		return wkF0 == 0 && wkF4 == 0 ? kdefault : (kaF0 * wkF0 + kaF4 * wkF4) / (wkF0 + wkF4);
	}

	/* Correct for multiple substitutions */
	CorrectKappaK80(P, Q) {
		let a = 1 - 2 * P - Q, b = 1 - 2 * Q;
		a = -Math.log(a); b = -Math.log(b);
		return 2 * a / b - 1;
	}
	CorrectKappaF84(P, Q, pi4) {
		let maxkappa = 2;

		let Y = pi4[0] + pi4[1], R = pi4[2] + pi4[3],
			tc = pi4[0] * pi4[1], ag = pi4[2] * pi4[3];

		//HKY85
		let A = tc / Y + ag / R, B = tc + ag, C = Y * R,
			a = B / A + (2 * (tc * R / Y + ag * Y / R) * (1 - Q / (2 * C)) - P) / (2 * A),
			b = 1 - Q / (2 * C);
		a = -Math.log(a);
		b = -Math.log(b);
		let kappa = 1 + A / B * (a / b - 1)
		if (kappa > 99) kappa = maxkappa
		return kappa; /* k_F84=>k_HKY85 */
	}
	CorrectKaksJC69(N, P, Q) {
		let r = (P + Q >= 0.75) ? 0.75 * (N - 1) / N : P + Q
		return -0.75 * Math.log(1 - r * 4 / 3);
	}
	CorrectKaksK80(P, Q) {
		let a = 1 - 2 * P - Q, b = 1 - 2 * Q;
		a = -Math.log(a); b = -Math.log(b);
		return 0.5 * a + 0.25 * b;
	}
	CorrectKaksF84(P, Q, pi4) {
		let Y = pi4[0] + pi4[1], R = pi4[2] + pi4[3],
			tc = pi4[0] * pi4[1], ag = pi4[2] * pi4[3];

		//HKY85
		let A = tc / Y + ag / R, B = tc + ag, C = Y * R,
			a = 1 - Y * R * P / (2 * (tc * R + ag * Y)) - Q / 2 * (tc * R / Y + ag * Y / R) / (tc * R + ag * Y),
			b = 1 - Q / (2 * Y * R);

		a = - Math.log(a);
		b = - Math.log(b);

		return 2 * (ag * a / R + tc * a / Y + (C - tc * R / Y - ag * Y / R) * b);
	}

	/**
	* Count synonymous(Sd) and nonsynonymous(Nd) differences
	* @param {string} seq1 encoded DNA sequence
	* @param {string} seq2 encoded DNA sequence
	* @param {string} pi codon frequency
	* @param {string} kappa estimated substitution rate
	* @param {string} omega estimated Ka/Ks
	* @param {string} t estimated divergence time
	* @returns {[number, number, number, number]} `[Sdts, Sdtv, Ndts, Ndtv]` number of synonymous and nonsynonymous substitutions in transition and transversion sites
	*/
	CountDiffs(seq1, seq2, t, omega, kappa, pi) {
		let Sdts = 0, Sdtv = 0,
			Ndts = 0, Ndtv = 0

		let depi = [], k = 0;
		pi.forEach((x, i) => { depi[i] = k; if (x > 0) k++ })
		let enpi = [];
		pi.forEach((x, i) => { if (x > 0) enpi[enpi.length] = i })
		pi = pi.filter(x => x > 0)

		let PMatrix = this.GetPMatCodon(t, omega, kappa, pi, enpi);

		let arr1 = seq1.match(/.{3}/g),
			arr2 = seq2.match(/.{3}/g);
		for (let h in arr1) {
			let codon1 = arr1[h],
				codon2 = arr2[h];
			if (codon1 == codon2) continue;

			/* syn ts & tv, nonsyn ts & tv for 2 codons */
			let sts = 0, stv = 0, nts = 0, ntv = 0;
			//diff[]: position of different codon 
			let diff = [0, 1, 2].filter(i => codon1[i] != codon2[i]),
				//ndiff: differences of two codons
				ndiff = diff.length,
				npath = { 1: 1, 2: 2, 3: 6 }[ndiff];

			if (ndiff == 1) {
				let transi = this.isTransition(codon1[diff[0]], codon2[diff[0]]);
				if (this.isSynonyms(codon1, codon2)) {
					if (transi) sts++; else stv++;
				} else {
					if (transi) nts++; else ntv++;
				}
			}
			else {   /* ndiff=2 or 3 */
				let nstop = 0,
					ppath = new Array(npath).fill(1),
					stspath = new Array(npath).fill(0),
					stvpath = new Array(npath).fill(0),
					ntspath = new Array(npath).fill(0),
					ntvpath = new Array(npath).fill(0);
				for (let k = 0; k < npath; k++) {
					//set the step[] 
					let step = [];
					if (ndiff == 2) {
						step[0] = diff[k];
						step[1] = diff[1 - k];
					}
					else {
						step[0] = k / 2 | 0;
						step[1] = k % 2;
						if (step[0] <= step[1])
							step[1]++;
						step[2] = 3 - step[0] - step[1];
					}//end of set the step[]

					// ct1, ct2: mutated codon
					let ct1 = codon1, ct2 = codon1;

					//ppath[]: probabilty of each path
					for (let i1 = 0; i1 < ndiff; i1++) {
						ct2 = this.mutate(ct2, step[i1], codon2[step[i1]])
						//ppath[k]: probabilty of path k
						ppath[k] *= PMatrix[depi[this.encodeCodon(ct1)]][depi[this.encodeCodon(ct2)]];

						if (this.getAminoAcid(ct1) == '!') {
							nstop++;
							ppath[k] = 0;
							break;
						}

						let transi = this.isTransition(codon1[step[i1]], codon2[step[i1]])

						//ts & tr when syn & nonsyn in path k
						if (this.isSynonyms(ct1, ct2)) {
							if (transi) stspath[k]++;
							else stvpath[k]++;
						}
						else {
							if (transi) ntspath[k]++;
							else ntvpath[k]++;
						}
						ct1 = ct2;
					}

				}  /* for(k,npath) */
				if (npath == nstop) {  /* all paths through stop codons */
					if (ndiff == 2) {
						nts = 0.5; ntv = 1.5;
					}
					else {
						nts = 0.5; ntv = 2.5;
					}
				}
				else {
					//sum probabilty of all path
					let sump = VectorUtils.sum(ppath);
					if (sump > 1e-20) {
						for (let k = 0; k < npath; k++) { //p: the probabilty of path k
							let p = ppath[k] / sump;
							sts += stspath[k] * p; stv += stvpath[k] * p;
							nts += ntspath[k] * p; ntv += ntvpath[k] * p;
						}
					}
				}
			}//end of if(ndiff)
			Sdts += sts; Sdtv += stv;
			Ndts += nts; Ndtv += ntv;
		}//end of for(h)

		return [Sdts, Sdtv, Ndts, Ndtv];
	}

	/**
	* Calculate transition probability matrix
	* @param {string} t estimated divergence time
	* @param {string} omega estimated Ka/Ks
	* @param {string} kappa estimated substitution rate
	* @param {string} pi codon frequency
	* @param {string} enpi codon frequency convertion array
	* @returns {number[][]} transition probability matrix
	*/
	GetPMatCodon(t, omega, kappa, pi, enpi) {
		let PMatrix = new Array(pi.length).fill(0).map(_ => []);
		for (let i = 0; i < pi.length; i++) {
			for (let j = 0; j < i; j++) {
				//codon 'from'
				let codon1 = this.decodeCodon(enpi[i]);
				//codon 'to'
				let codon2 = this.decodeCodon(enpi[j]);
				//stop codon
				if (this.getAminoAcid(codon1) == '!' ||
					this.getAminoAcid(codon2) == '!') {
					PMatrix[j][i] = PMatrix[i][j] = 0;
					continue;
				}

				//whether two codons only have one difference
				let diff = [0, 1, 2].filter(i => codon1[i] != codon2[i])
				if (diff.length == 1) {
					//only have one difference
					let r = 1;
					//transition
					if (this.isTransition(codon1[diff[0]], codon2[diff[0]]))
						r *= kappa;

					//nonsynonymous
					if (!this.isSynonyms(codon1, codon2))
						r *= omega;

					//diagonal element is equal
					PMatrix[j][i] = PMatrix[i][j] = r;
				} else {
					PMatrix[j][i] = PMatrix[i][j] = 0
				}
			}
		}

		for (let i in PMatrix) {
			//PMatrix[](*Q): transition probability matrix
			PMatrix[i] = PMatrix[i].map((x, j) => x * pi[j]);
			//scale the sum of PMat[][j](j=0-63) to zero
			PMatrix[i][i] = -VectorUtils.sum(PMatrix[i]);
		}
		//The sum of transition probability of main diagnoal elements
		let mr = pi.reduce((a, p, i) => a - p * PMatrix[i][i], 0)

		//calculate exp(PMatrix*t)
		let [Root, U, V] = this.eigenQREV(PMatrix, pi);
		for (let i in Root) Root[i] /= mr;
		PMatrix = this.PMatUVRoot(t, U, V, Root);
		return PMatrix;
	}

	/* P(t) = U * exp{Root*t} * V */
	PMatUVRoot(t, U, V, Root) {
		let n = Root.length;
		let smallp = 0;
		let P = new Array(n).fill(0).map(_ => new Array(n).fill(0))

		for (let k = 0; k < n; k++)
			for (let i = 0, expt = Math.exp(t * Root[k]); i < n; i++)
				for (let j = 0, uexpt = U[i][k] * expt; j < n; j++)
					P[i][j] += uexpt * V[k][j];

		P = P.map(x => x.map(x => x < smallp ? 0 : x))

		return P;
	}

	/**
	* Count synonymous(S) and nonsynonymous(N) sites
	* @param {string} seq encoded DNA sequence
	* @param {string} pi codon frequency
	* @param {string} kappa estimated substitution rate
	* @returns {[number, number]} `[S, N]` number of synonymous and nonsynonymous sites
	*/
	CountSites(seq, pi, kappa) {
		let N = 0, S = 0,
			fbS = [0, 0, 0, 0], //probablity of syn of four nul.
			fbN = [0, 0, 0, 0] //probablity of non-syn of four nul.

		for (let codon of seq.match(/.{3}/g)) {
			for (let pos of [0, 1, 2]) {
				for (let c of this.getMutatedCodonList(codon, pos)) {
					if (this.getAminoAcid(c) == "!") continue;
					let r = pi[this.encodeCodon(c)];
					if (this.isTransition(c[pos], codon[pos]))
						r *= kappa;
					if (this.isSynonyms(codon, c)) {
						fbS[codon[pos]] += r
					} else {
						fbN[codon[pos]] += r
					}
				}
			}
		};

		S = VectorUtils.sum(fbS);
		fbS = fbS.map(x => x / S);
		N = VectorUtils.sum(fbN);
		fbN = fbN.map(x => x / N);
		//Scale Stot+Ntot to seq.length
		let r = seq.length / (S + N);
		S *= r; N *= r;

		return [S, N, fbS, fbN];
	}

	/**
	* calculate Ka/Ks
	* @param {string} seq1 encoded DNA sequence1
	* @param {string} seq2 encoded DNA sequence2
	* @returns {{Ka: number, Ks: number, S: number, N: number, Sd: number, Nd: number, t: number, method: string, KAPPA: number[]}}
	*/
	calc(seq1, seq2) {
		let pi = this.getFreqency(seq1 + seq2);
		let kappa = this.GetKappa(seq1, seq2);

		let fbS = new Array(4).fill(0),
			fbN = new Array(4).fill(0);
		let S = 0, N = 0, Sd, Nd;

		//Count sites of sequence 1
		let [St, Nt, fbSt, fbNt] = this.CountSites(seq1, pi, kappa);
		S += St / 2; N += Nt / 2;
		fbSt.forEach((x, i) => fbS[i] += x / 2);
		fbNt.forEach((x, i) => fbN[i] += x / 2);

		//Count sites of sequence 2
		[St, Nt, fbSt, fbNt] = this.CountSites(seq2, pi, kappa);
		S += St / 2; N += Nt / 2;
		fbSt.forEach((x, i) => fbS[i] += x / 2);
		fbNt.forEach((x, i) => fbN[i] += x / 2);

		let nround = 10, dS, dN;
		let w0 = 0, dS0 = 0, dN0 = 0, omega = 1, t = 0.4,
			accu = 5e-4, minomega = 1e-5, maxomega = 99;
		for (let ir = 0; ir < nround; ir++) {   /* iteration */
			let [Sdts, Sdtv, Ndts, Ndtv] = this.CountDiffs(seq1, seq2, t, omega, kappa, pi);

			Sd = Sdts + Sdtv;
			Nd = Ndts + Ndtv;

			//synonymous
			dS = this.CorrectKaksF84(Sdts / S, Sdtv / S, fbS) || this.CorrectKaksK80(Sdts / S, Sdtv / S) || this.CorrectKaksJC69(S, Sdts / S, Sdtv / S);
			//nonsynonymous
			dN = this.CorrectKaksF84(Ndts / N, Ndtv / N, fbN) || this.CorrectKaksK80(Ndts / N, Ndtv / N) || this.CorrectKaksJC69(N, Ndts / N, Ndtv / N);

			if (Math.abs(dS - dS0) < accu && Math.abs(dN - dN0) < accu && Math.abs(omega - w0) < accu)
				break;

			t = dS * 3 * S / (S + N) + dN * 3 * N / (S + N);
			omega = dS < 1e-9 ? maxomega : Math.max(minomega, dN / dS);

			[dS0, dN0, w0] = [dS, dN, omega];
		} //end of for(ir) */

		return {
			Ks: dS, Ka: dN, S, N, Sd, Nd,
			t: (S * dS + N * dN) / (S + N),
			KAPPA: [kappa, kappa, 1, 1, 1, 1],
			method: this.method
		}
	}

	//The following functions are used to calculate PMatrix by Taylor.

	/**
	* Finds the eigen solution of the rate matrix Q for a time-reversible Markov process, using the algorithm for a real symmetric matrix.
	*
	* Ziheng Yang, 25 December 2001
	* @param {string} Q symmetric transition probability matrix
	* @param {string} pi codon frequency
	* @returns {[number[][], number[], number[]]} `[Root, U, V]`, where `[U] [Q] [U^-1] = [Root]`
	*/
	eigenQREV(Q, pi) {
		let pi_sqrt = pi.map(x => Math.sqrt(x));

		let A = Q.map((v, i) => v.map((x, j) => x * pi_sqrt[i] / pi_sqrt[j]));

		let [U, Root] = this.eigenRealSym(A);
		let V = new Array(U.length).fill(0).map(_ => []);
		for (let i in U)
			for (let j in U)
				V[i][j] = U[j][i] * pi_sqrt[j];
		U = U.map((v, i) => v.map(x => x / pi_sqrt[i]));
		return [Root, U, V];
	}

	HouseholderRealSym(a) {
		let n = a.length
		let d = [], e = []

		for (let i = n - 1; i >= 1; i--) {
			let m = i - 1;
			let h = 0, scale = 0;
			if (m > 0) {
				scale = VectorUtils.sum(a[i].slice(0, m));
				if (scale == 0)
					e[i] = a[i][m];
				else {
					h = a[i].slice(0, m + 1).reduce((a, x) => a + (x / scale) ** 2, 0);
					for (let k = 0; k <= m; k++)
						a[i][k] /= scale;
					let t = a[i][m];
					let g = (t >= 0 ? -Math.sqrt(h) : Math.sqrt(h));
					e[i] = scale * g;
					h -= t * g;
					a[i][m] = t - g;
					let f = 0;
					for (let j = 0; j <= m; j++) {
						a[j][i] = a[i][j] / h;
						let g = 0;
						for (let k = 0; k <= j; k++)
							g += a[j][k] * a[i][k];
						for (let k = j + 1; k <= m; k++)
							g += a[k][j] * a[i][k];
						e[j] = g / h;
						f += e[j] * a[i][j];
					}
					let hh = f / (h * 2);
					for (let j = 0; j <= m; j++) {
						let f = a[i][j];
						e[j] = g = e[j] - hh * f;
						for (let k = 0; k <= j; k++)
							a[j][k] -= (f * e[k] + g * a[i][k]);
					}
				}
			}
			else e[i] = a[i][m];
			d[i] = h;
		}
		d[0] = e[0] = 0;

		/* Get eigenvectors */
		for (let i = 0; i < n; i++) {
			let m = i - 1;
			if (d[i]) {
				for (let j = 0; j <= m; j++) {
					let g = 0;
					for (let k = 0; k <= m; k++)
						g += a[i][k] * a[k][j];
					for (let k = 0; k <= m; k++)
						a[k][j] -= g * a[k][i];
				}
			}
			d[i] = a[i][i];
			a[i][i] = 1;
			for (let j = 0; j <= m; j++) a[j][i] = a[i][j] = 0
		}

		return [d, a, e]
	}

	/**
	* Finds the eigen solution of a tridiagonal matrix represented by eigen values `d` and off-diagonal `e`.
	*
	* Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK. Ziheng Yang, 23 May 2001.
	* @param {number[]} d eigen values
	* @param {number[]} e off-diagonal matrix
	* @param {number[][]} U the eigen solution to the original real symmetric matrix
	* @returns {[number[], number[][]]} `[d, U]`, right vectors and eigenvalues
	*/
	EigenTridagQLImplicit(d, e, z) {
		let m, iter, niter = 30, i, k;
		let s, r, p, g, f, dd, c, b, aa, bb, n = d.length;

		for (let i = 1; i < n; i++) e[i - 1] = e[i]; e[n - 1] = 0;
		for (let j = 0; j < n; j++) {
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

					if ((aa = Math.abs(g)) > 1) r = aa * Math.sqrt(1 + 1 / (g * g));
					else r = Math.sqrt(1 + g * g);

					g = d[m] - d[j] + e[j] / (g + (g >= 0 ? 1 : -1) * Math.abs(r));
					s = c = 1;
					p = 0;
					for (i = m - 1; i >= j; i--) {
						f = s * e[i];
						b = c * e[i];

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
							f = z[k][i + 1];
							z[k][i + 1] = s * z[k][i] + c * f;
							z[k][i] = c * z[k][i] - s * f;
						}
					}
					if (r == 0 && i >= j) continue;
					d[j] -= p; e[j] = g; e[m] = 0;
				}
			} while (m != j);
		}
		return [d, z];
	}

	/**
	* Sorts the eigen values `d` and rearrange the right eigen vectors `U`.
	* @param {number[]} d eigen values
	* @param {number[][]} U right eigen vectors
	* @returns {[number[], number[][]]} `[d, U]`, sorted right vectors and eigenvalues
	*/
	EigenSort(d, U) {
		let k, j, i;
		let p, n = d.length;

		for (i = 0; i < n - 1; i++) {
			p = d[k = i];
			for (j = i + 1; j < n; j++)
				if (d[j] >= p) p = d[k = j];
			if (k != i) {
				d[k] = d[i];
				d[i] = p;
				for (j = 0; j < n; j++) {
					p = U[j][i];
					U[j][i] = U[j][k];
					U[j][k] = p;
				}
			}
		}
		return [d, U]
	}

	/**
	* Finds the eigen solution of a real symmetrical matrix A[n*n].
	* The matrix is first reduced to a tridiagonal matrix using `HouseholderRealSym()`,
	* and then using the QL algorithm with implicit shifts.
	*
	* Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK. Ziheng Yang, 23 May 2001.
	* @param {number[][]} A symmetrical matrix
	* @returns {[number[][], number[]]} `[A, Root]`, right vectors and eigenvalues
	*/
	eigenRealSym(A) {
		let Root, temp;
		[Root, A, temp] = this.HouseholderRealSym(A);
		[Root, A] = this.EigenTridagQLImplicit(Root, temp, A);
		[Root, A] = this.EigenSort(Root, A);
		return [A, Root];
	}
}

class MYN extends YN {
	constructor(genetic_code = 1) {
		super(genetic_code);
		this.method = "MYN";
	}

	/**
	* Estimate substitution rate (kappa) using the fourfold degenerate sites at third codon positions and nondegenerate sites
	* @param {string} seq1 encoded DNA sequence
	* @param {string} seq2 encoded DNA sequence
	* @returns {number, number} `[kappatc, kappaag]` T-C and A-G substitution rate
	*/
	GetKappa(seq1, seq2) {
		let kdefault = 2;
		let arr1 = seq1.match(/.{3}/g),
			arr2 = seq2.match(/.{3}/g);

		let F0 = new Array(4).fill(0).map(_ => new Array(4).fill(0)),
			F4 = new Array(4).fill(0).map(_ => new Array(4).fill(0));
		//Get Pi[] of A,C,G,T

		for (let h in arr1) {
			let codon1 = arr1[h],
				codon2 = arr2[h];

			//Find non-degenerate sites
			for (let pos of [0, 1, 2])
				if (
					this.getMutatedCodonList(codon1, pos)
						.every(c => !this.isSynonyms(c, codon1)) &&
					this.getMutatedCodonList(codon2, pos)
						.every(c => !this.isSynonyms(c, codon2))
				) {
					F0[codon1[pos]][codon2[pos]] += .5;
					F0[codon2[pos]][codon1[pos]] += .5;
				}

			// Find 4-fold degenerate sites at 3rd position
			if (
				this.isSynonyms(codon1, codon2) &&
				this.getMutatedCodonList(codon1, 2)
					.every(c => this.isSynonyms(c, codon1)) &&
				this.getMutatedCodonList(codon2, 2)
					.every(c => this.isSynonyms(c, codon2))
			) {
				F4[codon1[2]][codon2[2]] += .5;
				F4[codon2[2]][codon1[2]] += .5;
			}
		}//end of for(h)	

		let freq2Kappa = (F) => {
			let S = VectorUtils.sum(F.flat());
			if (S < 0) return 0;
			F = F.map(x => x.map(y => y / S))
			//Transition
			let T1 = F[2][3] * 2, T2 = F[0][1] * 2;
			//Tranversion
			let V = 1 - T1 - T2 - [0, 1, 2, 3].reduce((a, i) => a + F[i][i], 0);
			//pi[]: the sum probabilty of T, C, A, G, respectively
			let pi4 = [0, 1, 2, 3].map(i => VectorUtils.sum(F[i]));
			//Correct kappa
			let [katc, kaag] = this.CorrectKappaTN93(T1, T2, V, pi4);
			let wk = (katc > 0 && kaag > 0) ? S : 0
			return [katc, kaag, wk];
		}

		let [katcF0, kaagF0, wkF0] = freq2Kappa(F0),
			[katcF4, kaagF4, wkF4] = freq2Kappa(F4);

		if (wkF0 == 0 && wkF4 == 0)
			return [kdefault, kdefault]
		else {
			let kappatc = (katcF0 * wkF0 + katcF4 * wkF4) / (wkF0 + wkF4),
				kappaag = (kaagF0 * wkF0 + kaagF4 * wkF4) / (wkF0 + wkF4)
			return [kappatc, kappaag]
		}
	}

	/* Correct kappas */
	CorrectKappaTN93(P1, P2, Q, pi4) {
		let Y = pi4[0] + pi4[1], R = pi4[2] + pi4[3],
			tc = pi4[0] * pi4[1], ag = pi4[2] * pi4[3];

		let A = tc / Y + ag / R, B = tc + ag, C = Y * R,
			a1 = 1 - R * P1 / (2 * ag) - Q / (2 * R),
			a2 = 1 - Y * P2 / (2 * tc) - Q / (2 * Y),
			b = 1 - Q / (2 * C);

		a1 = -Math.log(a1);
		a2 = -Math.log(a2);
		b = -Math.log(b);
		//Kappa
		return [(a2 / b - R) / Y, (a1 / b - Y) / R]
	}

	/* Correct Ka and Ks */
	CorrectKaksTN93(P1, P2, Q, pi4) {
		let Y = pi4[0] + pi4[1], R = pi4[2] + pi4[3],
			tc = pi4[0] * pi4[1], ag = pi4[2] * pi4[3];

		let A = tc / Y + ag / R, B = tc + ag, C = Y * R,
			a1 = 1 - R * P1 / (2 * ag) - Q / (2 * R),
			a2 = 1 - Y * P2 / (2 * tc) - Q / (2 * Y),
			b = 1 - Q / (2 * Y * R);

		a1 = -Math.log(a1);
		a2 = -Math.log(a2);
		b = -Math.log(b);
		//Ka or Ks
		return 2 * ((ag * a1 / R) + (tc * a2 / Y) + ((C - ag * Y / R - tc * R / Y) * b));
	}

	/**
	* Calculate transition probability matrix
	* @param {string} t estimated divergence time
	* @param {string} omega estimated Ka/Ks
	* @param {string} kappa estimated substitution rate
	* @param {string} pi codon frequency
	* @param {string} enpi codon frequency convertion array
	* @returns {number[][]} transition probability matrix
	*/
	GetPMatCodon(t, omega, kappatc, kappaag, pi, enpi) {
		let PMatrix = new Array(pi.length).fill(0).map(_ => []);
		for (let i = 0; i < pi.length; i++) {
			for (let j = 0; j < i; j++) {
				//codon 'from'
				let codon1 = this.decodeCodon(enpi[i]);
				//codon 'to'
				let codon2 = this.decodeCodon(enpi[j]);
				//stop codon
				if (this.getAminoAcid(codon1) == '!' ||
					this.getAminoAcid(codon2) == '!') {
					PMatrix[j][i] = PMatrix[i][j] = 0;
					continue;
				}

				let diff = [0, 1, 2].filter(i => codon1[i] != codon2[i])
				//whether two codons only have one difference
				if (diff.length == 1) {
					//only have one difference
					let r = 1;
					//transition
					let transi = +codon1[diff[0]] + +codon2[diff[0]]
					if (transi == 1)
						r *= kappatc;//T<->C
					else if (transi == 5)
						r *= kappaag;//A<->G
					//nonsynonymous
					if (!this.isSynonyms(codon1, codon2))
						r *= omega;

					//diagonal element is equal
					PMatrix[j][i] = PMatrix[i][j] = r;
				} else {
					PMatrix[j][i] = PMatrix[i][j] = 0
				}
			}
		}

		for (let i in PMatrix) {
			//PMatrix[](*Q): transition probability matrix
			PMatrix[i] = PMatrix[i].map((x, j) => x * pi[j]);
			//scale the sum of PMat[][j](j=0-63) to zero
			PMatrix[i][i] = -VectorUtils.sum(PMatrix[i]);
		}
		//The sum of transition probability of main diagnoal elements
		let mr = pi.reduce((a, p, i) => a - p * PMatrix[i][i], 0)

		//calculate exp(PMatrix*t)
		let [Root, U, V] = this.eigenQREV(PMatrix, pi);
		for (let i in Root) Root[i] /= mr;
		PMatrix = this.PMatUVRoot(t, U, V, Root);
		return PMatrix;
	}

	/**
	* Count synonymous(Sd) and nonsynonymous(Nd) differences
	* @param {string} seq1 encoded DNA sequence
	* @param {string} seq2 encoded DNA sequence
	* @param {string} pi codon frequency
	* @param {string} kappatc estimated T-C substitution rate
	* @param {string} kappaag estimated A-G substitution rate
	* @param {string} omega estimated Ka/Ks
	* @param {string} t estimated divergence time
	* @returns {[number, number, number, number, number, number]} `[Sdts1, Sdts2, Sdtv, Ndts1, Ndts2, Ndtv]` number of synonymous and nonsynonymous substitutions in transition and transversion sites
	*/
	CountDiffs(seq1, seq2, t, omega, kappatc, kappaag, pi) {
		let Sdts1 = 0, Sdts2 = 0, Sdtv = 0,
			Ndts1 = 0, Ndts2 = 0, Ndtv = 0;

		let depi = [], k = 0;
		pi.forEach((x, i) => { depi[i] = k; if (x > 0) k++ })
		let enpi = [];
		pi.forEach((x, i) => { if (x > 0) enpi[enpi.length] = i })
		pi = pi.filter(x => x > 0)

		let PMatrix = this.GetPMatCodon(t, omega, kappatc, kappaag, pi, enpi);

		let arr1 = seq1.match(/.{3}/g),
			arr2 = seq2.match(/.{3}/g);
		for (let h in arr1) {
			let codon1 = arr1[h],
				codon2 = arr2[h];
			if (codon1 == codon2) continue;
			/* syn ts & tv, nonsyn ts & tv for 2 codons */
			let sts1 = 0, sts2 = 0, stv = 0,
				nts1 = 0, nts2 = 0, ntv = 0;
			//diff[]: position of different codon 
			let diff = [0, 1, 2].filter(i => codon1[i] != codon2[i]),
				//ndiff: differences of two codons
				ndiff = diff.length,
				npath = { 1: 1, 2: 2, 3: 6 }[ndiff];

			if (ndiff == 1) {
				let transi = +codon1[diff[0]] + +codon2[diff[0]];
				//transi=(transi==1 || transi==5);
				if (this.isSynonyms(codon1, codon2)) {
					if (transi == 5) sts1++;
					else if (transi == 1) sts2++;
					else stv++;
				}
				else {
					if (transi == 5) nts1++;
					else if (transi == 1) nts2++;
					else ntv++;
				}
			}
			else {   /* ndiff=2 or 3 */
				let nstop = 0,
					ppath = new Array(npath).fill(1),
					sts1path = new Array(npath).fill(0),
					sts2path = new Array(npath).fill(0),
					stvpath = new Array(npath).fill(0),
					nts1path = new Array(npath).fill(0),
					nts2path = new Array(npath).fill(0),
					ntvpath = new Array(npath).fill(0);
				for (let k = 0; k < npath; k++) {
					//set the step[]
					let step = []
					if (ndiff == 2) {
						step[0] = diff[k];
						step[1] = diff[1 - k];
					}
					else {
						step[0] = k / 2 | 0;
						step[1] = k % 2;
						if (step[0] <= step[1]) step[1]++;
						step[2] = 3 - step[0] - step[1];
					}//end of set the step[]

					let ct1 = codon1, ct2 = codon1;

					//ppath[]: probabilty of each path
					for (let i1 = 0; i1 < ndiff; i1++) {
						ct2 = this.mutate(ct2, step[i1], codon2[step[i1]])

						//ppath[k]: probabilty of path k
						ppath[k] *= PMatrix[depi[this.encodeCodon(ct1)]][depi[this.encodeCodon(ct2)]];

						if (this.getAminoAcid(ct1) == '!') {
							nstop++;
							ppath[k] = 0;
							break;
						}

						let transi = +codon1[step[i1]] + +codon2[step[i1]];

						//ts & tr when syn & nonsyn in path k
						if (this.isSynonyms(ct1, ct2)) {
							if (transi == 5) sts1path[k]++;
							else if (transi == 1) sts2path[k]++;
							else stvpath[k]++;
						}
						else {
							if (transi == 5) nts1path[k]++;
							else if (transi == 1) nts2path[k]++;
							else ntvpath[k]++;
						}
						ct1 = ct2;
					}
				}  /* for(k,npath) */
				if (npath == nstop) {  /* all paths through stop codons */
					if (ndiff == 2) {
						nts1 = 0.25; nts2 = 0.25; ntv = 1.5;
					}
					else {
						nts1 = 0.25; nts2 = 0.25; ntv = 2.5;
					}
				}
				else {
					//sum probabilty of all path
					let sump = VectorUtils.sum(ppath);
					if (sump > 1e-20) {
						for (let k = 0; k < npath; k++) { //p: the probabilty of path k
							let p = ppath[k] / sump;
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

	/**
	* Count synonymous(S) and nonsynonymous(N) sites
	* @param {string} seq encoded DNA sequence
	* @param {string} pi codon frequency
	* @param {string} kappatc estimated T-C substitution rate
	* @param {string} kappaag estimated A-G substitution rate
	* @returns {[number, number]} `[S, N]` number of synonymous and nonsynonymous sites
	*/
	CountSites(seq, pi, kappatc, kappaag) {
		let fbS = [0, 0, 0, 0], //probablity of syn of four nul.
			fbN = [0, 0, 0, 0] //probablity of non-syn of four nul.

		for (let codon of seq.match(/.{3}/g)) {
			for (let pos of [0, 1, 2]) {
				for (let c of this.getMutatedCodonList(codon, pos)) {
					if (this.getAminoAcid(c) == "!") continue;
					let r = pi[this.encodeCodon(c)];
					let transi = +codon[pos] + +c[pos]
					if (transi == 1) r *= kappatc;
					else if (transi == 5) r *= kappaag;
					if (this.isSynonyms(codon, c)) {
						fbS[codon[pos]] += r; //syn probability of A,C,G,T
					}
					else {
						fbN[codon[pos]] += r; //nonsyn probability of A,C,G,T
					}
				}
			}
		}

		let S = VectorUtils.sum(fbS);
		fbS = fbS.map(x => x / S);
		let N = VectorUtils.sum(fbN);
		fbN = fbN.map(x => x / N);
		//Scale Stot+Ntot to seq.length
		let r = seq.length / (S + N);
		S *= r; N *= r;

		return [S, N, fbS, fbN];
	}

	/**
	* calculate Ka/Ks
	* @param {string} seq1 encoded DNA sequence1
	* @param {string} seq2 encoded DNA sequence2
	* @returns {{Ka: number, Ks: number, S: number, N: number, Sd: number, Nd: number, t: number, method: string, KAPPA: number[]}}
	*/
	calc(seq1, seq2) {
		let pi = this.getFreqency(seq1 + seq2);
		let [kappatc, kappaag] = this.GetKappa(seq1, seq2);

		let fbS = new Array(4).fill(0),
			fbN = new Array(4).fill(0);
		let S = 0, N = 0, Sd, Nd;

		//Count sites of sequence 1		
		let [St, Nt, fbSt, fbNt] = this.CountSites(seq1, pi, kappatc, kappaag);
		S += St / 2; N += Nt / 2;
		fbSt.forEach((x, i) => fbS[i] += x / 2);
		fbNt.forEach((x, i) => fbN[i] += x / 2);

		//Count sites of sequence 2
		[St, Nt, fbSt, fbNt] = this.CountSites(seq2, pi, kappatc, kappaag);
		S += St / 2; N += Nt / 2;
		fbSt.forEach((x, i) => fbS[i] += x / 2);
		fbNt.forEach((x, i) => fbN[i] += x / 2);

		let nround = 100, dS, dN;
		let w0 = 0, dS0 = 0, dN0 = 0, omega = 0.5, t = 0.09,
			accu = 5e-8, minomega = 1e-5, maxomega = 99;
		//Iterative loop
		for (let ir = 0; ir < nround; ir++) {   /* iteration */
			//Count differences
			let [Sdts1, Sdts2, Sdtv, Ndts1, Ndts2, Ndtv] = this.CountDiffs(seq1, seq2, t, omega, kappatc, kappaag, pi);

			//Synonymous(Sd) and nonsynonymous(Nd) differences
			Sd = Sdts1 + Sdts2 + Sdtv;
			Nd = Ndts1 + Ndts2 + Ndtv;

			//Seldom happen
			if (Sd > S) {
				Sdts1 *= S / Sd; Sdts2 *= S / Sd; Sdtv *= S / Sd;
			}
			if (Nd > N) {
				Ndts1 *= N / Nd; Ndts2 *= N / Nd; Ndtv *= N / Nd;
			}

			//Ks
			dS = this.CorrectKaksTN93(Sdts1 / S, Sdts2 / S, Sdtv / S, fbS) || this.CorrectKaksF84((Sdts1 + Sdts2) / S, Sdtv / S, fbS) || this.CorrectKaksK80((Sdts1 + Sdts2) / S, Sdtv / S) || this.CorrectKaksJC69(S, (Sdts1 + Sdts2) / S, Sdtv / S);
			//Ka
			dN = this.CorrectKaksTN93(Ndts1 / N, Ndts2 / N, Ndtv / N, fbN) || this.CorrectKaksF84((Ndts1 + Ndts2) / N, Ndtv / N, fbN) || this.CorrectKaksK80((Ndts1 + Ndts2) / N, Ndtv / N) || this.CorrectKaksJC69(N, (Ndts1 + Ndts2) / N, Ndtv / N);

			if (Math.abs(dS - dS0) < accu && Math.abs(dN - dN0) < accu && Math.abs(omega - w0) < accu)
				break;

			t = 3 * (dS * S + dN * N) / (S + N);
			omega = dS < 1e-9 ? maxomega : Math.max(minomega, dN / dS);

			[dS0, dN0, w0] = [dS, dN, omega];

		} //end of for(ir) */

		return {
			Ks: dS, Ka: dN, S, N, Sd, Nd,
			t: (S * dS + N * dN) / (S + N),
			KAPPA: [kappatc, kappaag, 1, 1, 1, 1],
			method: this.method
		}
	}
}
