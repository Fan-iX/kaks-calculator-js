# Kakscalculator.js

A javascript reimplementation of [KaKs_calculator](https://ngdc.cncb.ac.cn/biocode/tools/BT000001).

KaKs_Calculator is capable for calculating selective pressure on both coding and non-coding sequences. For coding sequences, it integrates several methods to calculate nonsynonymous (Ka) and synonymous (Ks) substitution rates; particularly, it adopts model selection and model averaging to include as many features as needed for accurately capturing evolutionary information in protein-coding sequences.

### Compatibility

 - [ES6 class](https://caniuse.com/es6-class) is required.

## Usage

```html
<script src="src/base.js"></script>
<script src="src/KaKs.js"></script>
<script src="src/NG.js"></script>
<script src="src/LWL.js"></script>
<script src="src/LPB.js"></script>
<script src="src/YN.js"></script>
<script src="src/GY.js"></script>
<script src="src/MSMA.js"></script>
```

```js
let code=
let kk = new KAKS(genetic_code) // code type, 1-33, see below
kk.Run(`NP_000006.1(1-57)
ATGGACATTGAAGCATATTTTGAAAGAATTGGCTATAAGAACTCTAGGAACAAATTG
ATGGACATCGAAGCATACTTTGAAAGGATTGGTTACAAGAACTCAGTGAATAAATTG`, // axt string
    ["MYN", "GY"],
    { gy_models: ["HKY", "TIMEF"] }
)
```

<details>
<summary>codon table</summary>

```
   TTT           TGA  CTG            ATA     AAA  AGG             GGG
    |             |    |              |       |    |               |
 1: FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Standard
 2: FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS!!VVVVAAAADDEEGGGG // Vertebrate Mitochondrial
 3: FFLLSSSSYY!!CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Yeast Mitochondrial
 4: FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Mold Mitochondrial
 5: FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG // Invertebrate Mitochondrial
 6: FFLLSSSSYYQQCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Ciliate, Dasycladacean and Hexamita
 9: FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG // Echinoderm and Flatworm Mitochondrial
10: FFLLSSSSYY!!CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Euplotid Nuclear
11: FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Bacterial and Plant Plastid
12: FFLLSSSSYY!!CC!WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Alternative Yeast Nuclear
13: FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG // Ascidian Mitochondrial
14: FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG // Alternative Flatworm Mitochondrial
15: FFLLSSSSYY!QCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Blepharisma Nuclear
16: FFLLSSSSYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Chlorophycean Mitochondrial
21: FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG // Trematode Mitochondrial
22: FFLLSS!SYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Scenedesmus obliquus mitochondrial
23: FF!LSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Thraustochytrium Mitochondrial
24: FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG // Rhabdopleuridae Mitochondrial
25: FFLLSSSSYY!!CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Candidate Division SR1 and Gracilibacteria
26: FFLLSSSSYY!!CC!WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Pachysolen tannophilus Nuclear
27: FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Karyorelict Nuclear
28: FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Condylostoma Nuclear
29: FFLLSSSSYYYYCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Mesodinium Nuclear
30: FFLLSSSSYYEECC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Peritrich Nuclear
31: FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG // Blastocrithidia Nuclear
33: FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG // Cephalodiscidae Mitochondrial UAA-Tyr
```

</details>

## Notice

Calculation results may be slightly different from the C++ version due to implementation

## Citation

[KaKs_Calculator 3.0: Calculating Selective Pressure On Coding And Non-Coding Sequences](http://dx.doi.org/10.1016/j.gpb.2021.12.002)

Zhang Z, 2022 Jun - Genomics, Proteomics & Bioinformatics
