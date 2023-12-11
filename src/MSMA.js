/************************************************************
* Filename: MSMA.js
* Abstract: Definition of model-selected and model-averged methods' classes.

* Version: js-2.0
* Author: Fanix
* Date: December.11, 2022

* Adapted from: KaKs_Calculator/MSMA.cpp and KaKs_Calculator/MSMA.h
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.
 
* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009

* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Apr. 2006

  References:
  Goldman N, Yang Z  (1994)  A codon-based model of nucleotide
  substitution for protein-coding DNA sequences. Mol. Biol. 
  Evol. 11:725-736.

  Posada, D. and Buckley, T.R. (2004) Model Selection and Model Averaging
  in Phylogenetics: Advantages of Akaike Information Criterion and Bayesian
  Approaches over Likelihood Ratio Tests, Syst. Biol., 53, 793-808.

  Sullivan, J. and Joyce, P. (2005) Model Selection in Phylogenetics, 
  Annual Review of Ecology, Evolution, and Systematics, 36, 445-466.  
*************************************************************/

class MS extends Base {
	constructor() {
		super();
		this.name = "MS";
	}

	/* Calculate Ka and Ks based on a given model, similar to the method of GY */
	selectModel(seq1, seq2, c_model, result4MA) {
		let tmp = {
			result: "",
			AICc: 0,
			freq: new Array(CODON),
			rate: new Array(NUMBER_OF_RATES),
			w: 0,
			t: 0
		};
		let zz = new GY94(c_model);

		tmp.result = zz.Run(seq1, seq2);
		tmp.AICc = zz.AICc;
		tmp.freq = zz.com.pi
		tmp.rate = zz.KAPPA
		tmp.w = zz.com.omega;
		tmp.t = 3. * zz.t;

		result4MA.push(tmp);
	}

	/* Choose the estimates under a model with smallest AICc */
	Run(seq1, seq2, result4MA, details) {

		let i, j, pos;
		let candidate_models = ["JC", "F81", "K2P", "HKY", "TNEF", "TN", "K3P", "K3PUF", "TIMEF", "TIM", "TVMEF", "TVM", "SYM", "GTR"];

		//Calculate Ka and Ks using 14 models
		for (i = 0; i < MODELCOUNT; i++)
		 this.selectModel(seq1, seq2, candidate_models[i], result4MA);

		//Choose the results under a model with smallest AICc
		for (pos = i = 0; i < result4MA.length; i++) {
			if (result4MA[i].AICc < result4MA[pos].AICc) pos = i;
		}

		//Calculate the AICc difference, substract the smallest AICc
		let diff = result4MA.map(x => x.AICc - result4MA[pos].AICc);

		//Compute Akaike weights
		let w = new Array(MODELCOUNT).fill(0), sum;
		for (sum = i = 0; i < MODELCOUNT; i++) {
			//akaike wights of each model
			for (j = 0; j < MODELCOUNT; j++) {
				let power = -0.5 * diff[j] - (-0.5 * diff[i]);
				//Avoid overflow
				if (power > 709) power = 700;
				else if (power < -709) power = -700;
				//Normal
				w[i] += Math.exp(power);
			}
			w[i] = 1. / w[i];
		}

		//Add Akaike weights to results
		let tmp = "";
		for (i = 0; i < MODELCOUNT; i++) {

			//Add akaike weights
			tmp = result4MA[i].result;
			j = tmp.lastIndexOf('\t');
			result4MA[i].result = tmp.substr(j, tmp.length - j);
			tmp = tmp.slice(0, j);
			j = tmp.lastIndexOf('\t');
			tmp = tmp.slice(0, j+1) + w[i];
			result4MA[i].result = tmp + result4MA[i].result;

			//Details on model selection
			details += result4MA[i].result;
		}

		//Results at "pos" is more reliable, replace 'method name' by "MS".
		tmp = result4MA[pos].result;
		i = tmp.indexOf('\t');
		j = tmp.indexOf('\t', i + 1);
		result4MA[pos].result = tmp.substr(0, i + 1) + this.name;
		result4MA[pos].result += tmp.substr(j, tmp.length - j);

		return result4MA[pos].result;
	}
}


class MA extends GY94 {
	constructor() {
		super();
		this.name = "MA";
		this.Small_Diff = 1e-6;
		this.w_rndu = 123456757;
		this.com.ns = 2;
		this.lnL = 0.0;

		this.com.icode = genetic_code - 1;
		if (this.com.icode > 11) this.com.icode = 0;
		this.com.ncode = this.Nsensecodon = 64 - this.getNumNonsense(this.com.icode);

		this.com.nkappa = 5;
		this.com.np = 2 + this.com.nkappa;
	}


	Run(seq1, seq2, result4MA) {
		let i, j, pos;

		//Find the smallest AICc
		for (pos = 0, i = 1; i < result4MA.length; i++) {
			if (result4MA[i].AICc < result4MA[pos].AICc) pos = i;
		}

		//Calculate the AICc difference, substract the smallest AICc
		let diff = [];
		for (i = 0; i < MODELCOUNT; i++) diff[i] = result4MA[i].AICc - result4MA[pos].AICc;

		//Compute Akake weights
		let w = new Array(MODELCOUNT).fill(0);
		for (i = 0; i < MODELCOUNT; i++) {
			//Avoid overflow
			for (j = 0; j < MODELCOUNT; j++) {
				let power = -0.5 * diff[j] - (-0.5 * diff[i]);
				if (power > 709) power = 700;
				else if (power < -709) power = -700;
				w[i] += Math.exp(power);
			}
			w[i] = 1. / w[i];
		}

		let I = new Array(MODELCOUNT).fill(0).map(_ => new Array(NUMBER_OF_RATES).fill(0));
		for (i = 0; i < MODELCOUNT; i++)
			for (j = 0; j < NUMBER_OF_RATES; j++)
				I[i][j] = 0;

		//JC, F81:	  rTC==rAG =rTA==rCG==rTG==rCA

		//K2P, HKY:   rTC==rAG!=rTA==rCG==rTG==rCA
		I[2][0] = 1; I[2][1] = 1;
		I[3][0] = 1; I[3][1] = 1;
		//TNEF,TN:    rTC!=rAG!=rTA==rCG==rTG==rCA
		I[4][0] = 1; I[4][1] = 1;
		I[5][0] = 1; I[5][1] = 1;
		//K3P, K3PUF: rTC==rAG!=rTA==rCG!=rTG==rCA
		I[6][0] = 1; I[6][1] = 1; I[6][2] = 1; I[6][3] = 1;
		I[7][0] = 1; I[7][1] = 1; I[7][2] = 1; I[7][3] = 1;
		//TIMEF, TIM: rTC!=rAG!=rTA==rCG!=rTG==rCA
		I[8][0] = 1; I[8][1] = 1; I[8][2] = 1; I[8][3] = 1;
		I[9][0] = 1; I[9][1] = 1; I[9][2] = 1; I[9][3] = 1;
		//TVMEF, TVM: rTC==rAG!=rTA!=rCG!=rTG!=rCA
		I[10][0] = 1; I[10][1] = 1; I[10][2] = 1; I[10][3] = 1; I[10][4] = 1;
		I[11][0] = 1; I[11][1] = 1; I[11][2] = 1; I[11][3] = 1; I[11][4] = 1;
		//SYM, GTR:   rTC!=rAG!=rTA!=rCG!=rTG!=rCA
		I[12][0] = 1; I[12][1] = 1; I[12][2] = 1; I[12][3] = 1; I[12][4] = 1;
		I[13][0] = 1; I[13][1] = 1; I[13][2] = 1; I[13][3] = 1; I[13][4] = 1;


		//Average parameters
		let para = new Array(8).fill(0);

		//Model-averaged time t
		for (j = 0; j < MODELCOUNT; j++) para[0] += (w[j] * result4MA[j].t);
		//para[0] /= sum;

		//Model-averaged Substitution Rates
		let sum = new Array(NUMBER_OF_RATES).fill(0);
		this.com.KAPPA = new Array(8).fill(0);

		for (i = 0; i < NUMBER_OF_RATES - 1; i++) {
			for (j = 0; j < MODELCOUNT; j++) {
				this.com.KAPPA[i] += (w[j] * I[j][i] * result4MA[j].rate[i]);
				sum[i] += (w[j] * I[j][i]);
			}
			if (sum[i] < 1e-50) this.com.KAPPA[i] = 1;
			else this.com.KAPPA[i] /= sum[i];
			para[i + 1] = this.com.KAPPA[i];
		}

		para[6] = 0;
		//Model-averaged omega w
		for (j = 0; j < MODELCOUNT; j++) para[6] += (w[j] * result4MA[j].w);
		//para[6] /= sum;

		//Model-averaged Codon Frequencies
		for (i = 0; i < CODON; i++) {
			this.com.pi[i] = 0.0;
			for (j = 0; j < MODELCOUNT; j++) this.com.pi[i] += (w[j] * result4MA[j].freq[i]);
			//com.pi[i] /= sum;
		}

		/* Preprocess */
		this.preProcess(seq1, seq2);

		//Calculate maximum likelihood score
		this.lnL = this.lfun2dSdN(para, this.com.np);

		//Copy subsitution rates
		this.KAPPA = this.com.KAPPA.map(x => x);

		//Compute Ka and Ks
		let PMat;
		[this.S, this.Ks, this.Ka, PMat] = this.EigenQc(1, para[0], this.KAPPA, this.com.omega);

		this.N = this.com.ls * 3 - this.S;
		this.Sd = this.Ks * this.S;
		this.Nd = this.Ka * this.N;
		let scale = (this.Sd + this.Nd) / this.snp;
		this.Sd /= scale;
		this.Nd /= scale;

		this.lnL = -this.lnL;

		this.t = para[0] / 3;

		//For the method of Model averaging, parameters' number is not specific.
		this.AICc = NA;

		return this.parseOutput();
	}


}
