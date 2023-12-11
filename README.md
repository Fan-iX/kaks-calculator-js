# Kakscalculator2.js

A javascript rewritten of [KaKs_calculator2](https://sourceforge.net/projects/kakscalculator2/).

KaKs_calculator2 is a toolkit of incorporating gamma series methods and sliding window strategies.

This branch is a nearly line-for-line JS rewritten of the C++ code form the github Fork [kakscalculator2](https://github.com/kullrich/kakscalculator2/).

### Compatibility

 - [ES6 class](https://caniuse.com/es6-class) is required.

## Usage

```html
<script src="src/base.js"></script>
<script src="src/KaKs.js"></script>
<script src="src/NG86.js"></script>
<script src="src/LWL85.js"></script>
<script src="src/LPB93.js"></script>
<script src="src/YN00.js"></script>
<script src="src/MYN.js"></script>
<script src="src/GY94.js"></script>
<script src="src/MSMA.js"></script>
```

```js
let kk = new KAKS()
kk.Run(`NP_000006.1(1-57)
ATGGACATTGAAGCATATTTTGAAAGAATTGGCTATAAGAACTCTAGGAACAAATTG
ATGGACATCGAAGCATACTTTGAAAGGATTGGTTACAAGAACTCAGTGAATAAATTG`, // axt string
    ["-C", 1, "-M", "GMYN"]
)
```

## Notice

Calculation results may be slightly different from the original C++ version due to floating-point precision (?)

## Citation

Please cite the selected following reference when using the KaKs_Calculator2.0:

Da-Peng Wang, Hao-Lei Wan, Song Zhang and Jun Yu. γ-MYN: a new algorithm for estimating Ka and Ks with consideration of variable substitution rates. Biology Direct 2009, 4:20.

Dapeng Wang, Song Zhang, Fuhong He, Jiang Zhu, Songnian Hu, Jun Yu. How Do Variable Substitution Rates Influence Ka and Ks Calculations?Genomics Proteomics Bioinformatics. 2009 Sep;7(3):116-27.

Dapeng Wang, Yubin Zhang, Zhang Zhang, Jiang Zhu, Jun Yu. KaKs_Calculator 2.0: a toolkit incorporating gamma-series methods and sliding window strategies.Genomics Proteomics Bioinformatics. 2010 Mar;8(1):77-80.
