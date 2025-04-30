# memesuite-lite

The tools in the [MEME suite](https://meme-suite.org/meme/) are foundational for many sequence-based analyses; MEME itself for discovering repeating patterns in sequences, FIMO for scanning these patterns against long sequences, and TOMTOM for scoring these patterns against each other. As we enter an age of large-scale and machine learning-based analyses, these tools can continue to be critical components because they answer fundamental questions but can be challenging to use in practice. Specifically, these tools require the use of intermediary text files (MEME-formatted for PWMs and FASTA formatted for sequences) followed by either using a command-line tool or a web portal.

memesuite-lite is a re-implementation of some of the algorithms in the MEME suite as Python functions that can be easily plugged into existing notebooks and code-bases. These implementations are well documented, multi-threaded, fast, and simple to use because they are implemented using numba. As an example of the speed, the TOMTOM implementation can be thousands of times faster than the MEME suite command-line tool due to severe inefficiencies in the code. 

The goal of memesuite-lite is twofold: (1) to scale these algorithms to millions of examples on modest consumer hardware, and (2) to make doing so extremely easy for a user. Accordingly, our implementations are meant to be scaled up and our focus is on finding new approximations or algorithmic tricks to accelerate these algorithms without losing too much precision (though losing some is not a problem).

### Installation

`pip install memelite`

### TOMTOM

[[tutorial](https://github.com/jmschrei/memesuite-lite/blob/main/tutorials/Tutorial_TOMTOM.ipynb)]

TOMTOM is an algorithm for calculating similarity scores between two sets of PWMs. Importantly, TOMTOM is not a similarity score itself, but rather the procedure for faithfully calculating p-values given a similarity score. Accordingly, it has two steps. First, consider all possible alignments between each query motif and the target motif and record the best score across all alignments. Second, convert this score to a p-value using the appropriate background distribution. This algorithm will be immensely useful in annotating the genome because these PWMs do not have to be frequency-based matrices like those we see in JASPAR. Rather, they can even be one-hot encoded sequences, allowing you to compare seqlets identified using a machine learning model against a database of known motifs and annotate them automatically. Or, they can be the attribution scores themselves at these seqlets being compared against a database of CWMs identified via a method like TF-MoDISco.

See the tutorial for more details, but here is an example of using the code to map a one-hot encoded sequence against the JASPAR motif database.

```python
from memelite import tomtom
from memelite.io import read_meme   # The same as tangermeme but returns numpy arrays

from tangermeme.utils import one_hot_encode

q = [one_hot_encode("ACGTGT").double()]

targets = read_meme("../../../common/JASPAR2024_CORE_non-redundant_pfms_meme.txt")
target_pwms = [pwm for pwm in targets.values()]

p, scores, offsets, overlaps, strands = tomtom(q, target_pwms)
```

By default, TOMTOM returns matrices with shape `n_queries x n_targets`. When these numbers are large, as in our goal of processing 1M queries by 1M targets in 1 minute, the full similarity matrix might not actually fit in memory. Fortunately, in these large-scale cases we usually do not actually care about the scores between all queries and all targets. Rather, we care about the scores between all queries and some number of the closest matches. Usually, this number can actually be quite small, such as in the case where we simply want to annotate each query with the closest (or closest five) targets.

If you want to only get these results for some number of nearest neighbors you can pass in the `n_neighbors` parameter and get a `n_queries x n_neighbors` matrix, saving a significant amount of memory.

```python
p, scores, offsets, overlaps, strands, idxs = tomtom(target_pwms, target_pwms, n_nearest=100)
```

Note that, in addition to the matrices of alignment statistics as before, this will also return a matrix of indexes showing which target indexes gave rise to these statistics. These are ordered such that the first index of `idxs` is the best hit and the second index is the second best hit, etc.

### FIMO

[[tutorial](https://github.com/jmschrei/memesuite-lite/blob/main/tutorials/Tutorial_FIMO.ipynb)]

FIMO is an algorithm for scanning a set of PWMs against one or more one-hot encoded sequences and finding the statistically significant hits. The algorithm is composed of two steps: scanning the PWMs just like one would apply a convolution, and converting these scores to p-values using a background distribution calculated exactly from the PWM. Although FIMO and TOMTOM both involve comparing a PWM against some other entity, they differ in the assumptions that they make about this second entity. TOMTOM assumes that the entity is short and that overhangs are not only likely but critically important to score correctly. FIMO assumes that the entity is long and so overhangs do not matter and, accordingly, does not care too much about scoring them correctly. See the TOMTOM tutorial for a more in-depth description of the two algorithms.

See the tutorial above for more details but here is an example of applying FIMO using a set of motifs written out in a MEME file against a set of sequences written out in a FASTA file.

```python
from memelite import fimo

hits = fimo("../tests/data/test.meme", "../tests/data/test.fa") 
```

In this case, `hits` will be a list of 12 pandas DataFrames with one DataFrame for each motif in the MEME file. 

Although we passed in filepaths above, the function signature is flexible. You can also pass in a dictionary for the motifs where the keys are motif names and the values are either numpy arrays or PyTorch tensors encoding the PWMs, and you can pass in a numpy array or PyTorch tensor for the sequences that will be scanned.
