!!!danger
    THIS IS A WORK IN PROGRESS

# Expression Tree Evaluation

This section documents the general algorithm used to take a set of inputs and evaluate them to produce a set of results given a user-defined set of operations. To describe the method more compactly, we introduce the following notation:

\(\alpha_i \in A\)
: An entity for which its parameters are known as a member of the set of all runtime inputs.

\(\overline\alpha_i \in \overline A\)
: The corresponding entity expressed in indeterminate form[^1].

\(\overline\beta_j = f_j(\overline\alpha_S, \overline\beta_T) \in \overline B, j \not\in T\)
: An intermediate result of the fuller computation as a function of some subset of \(A\) and \(B\).

\(\overline\gamma_k = g_k(\overline\alpha_S, \overline\beta_T) \in \overline\Gamma \)
: The result of a requested evaluation as a function of a subset of the input and intermediate expressions.

\(\overline\lambda: \overline A \rightarrow \overline\Gamma\)
: A user defined function which specifies how to map from a set of indeterminate inputs to a set of indeterminate results.

\(\lambda: A \rightarrow B \rightarrow \Gamma\)
: The computation's execution strategy as a map from the runtime-valued inputs to the runtime-valued outputs via a set of intermediate results. The map \(\lambda\) is induced by the user-defined map \(\overline\lambda\) which is specified literally as an anonymous lambda function (hence the notation). In words, a lambda defines what the user wishes to compute, and big-lambda defines how we intend on satisfying that wish.

As you can see, we adopt the convention that a bar overhead a symbol denotes a compile-time indeterminate expression, whereas the absense of the bar denotes a runtime-valued quantity. The goal is to find \(\lambda\) in a manner that provides strict upper bounds on memory usage and compiler runtime (within a factor of the complexity of the worst-case member of \(\Gamma\)). The challenges involved in identifying an optimal \(\lambda\) are numerous:

- Reification of the final results from \(\overline\alpha\) and \(\overline\beta\) to \(\gamma\) must happen in the final pass jointly or the compiler will fail to inline the code.
- Before expanding the expression tree, the complexity isn't fully known, and as such, the memory required to allocate space to accommodate \(B\) is similarly unknown.
- Processing an evaluation to produce an intermediate result generally increases as \(O(\prod_A |\overline\alpha|)\) (in words, the size of a resultant indeterminate expression scales with the product of the sizes of the inputs). Mismanaging this (i.e. going straight from \(A\) to \(\Gamma\)) can quickly explode and make compilation slow or impossible on finite hardware.
- By the same token, attempting too hard to identify an optimal \(\lambda\) at compile-time is in and of itself an expensive proposition.
- It is often the case that sub-expressions are dependencies of multiple results. In such a case, we would like to store the result as a \(\beta_j\). Even if the compiler can ultimately eliminate such an expression in a late-stage optimization pass, such elimination is not guaranteed depending on the complexity of \(\overline\lambda\) and also places an unnecessary burden on the compiler.
- Perhaps most-restrictive of all is that the user-defined computation supplied as an anonymous lambda function cannot be constructed or invoked as it is itself not treated by the compiler as a constant expression.[^2] This restriction implies that the expression trees cannot be described using flattened structs and arrays and must be encoded as types instead (i.e. the information the user can provide is based on the `std::invoke_result_t` of the input algorithm).

For the purposes of GAL, we relax the restriction that \(\lambda\) is necessarily optimal, and instead employ a number of heuristics to get a *close enough* result that is, at the very least, guaranteed to compile on most machines and produce accurate results. A high level description of the algorithm follows, split into several labeled phases:

### Expression tree construction

#### Input

The inputs are the set \(\alpha_i \in A\) and the user-defined function \(\overline\lambda\).

#### Algorithm

1. During expression tree construction, each operation node computes the checksum of its inputs and the operation. In the case of a unary operation, the checksum has two arguments (input checksum, and operation code), and in the case of a binary operation, the checksum has three arguments. In this manner, any two subtree checksums can be compared for equality in constant time. For leaf nodes, the checksum is fully determined by the id of the first indeterminate.[^3]
1. Each tree maintains at a top level, a list of additional trees representing the set \(\overline B\) of intermediate subexpressions identified.
1. Each node maintains a checksum combining itself and all of its children, making the tree a special form of a Merkle tree.[^4]
1. Additionally, the expression tree maintains a both continuously updated bloom filter.[^5] The advantages of the bloom filter is that it occupies constant memory which is extremely desirable for compile-time computation.
1. When a binary operation joins two nodes, it is possible that a child or the root of the LHS matches a child or a root of the RHS.
    1. Without loss of generality, swap the LHS and RHS if necessary so that the number of nodes in the LHS is less than the number of nodes in the RHS.
    1. To detect a common subexpression, traverse the LHS and perform an existence test for each node on the bloom filter of the RHS.[^6] Initial implementations may skip the initial bloom filtration and go straight to an \(O(n)\) scan of the tree (which may be sufficient for small inputs).[^7]
    1. If a match is found, we promote the subtree to the top level \(\overline B\), and replace both matching nodes with a placeholder node \(\overline\beta\). The node \(\overline\beta\) is given an identifier \(j\) corresponding to the index of the subtree in the top-level list.
    1. The top level lists of \(\overline B_\textit{LHS}\) and \(\overline B_\textit{RHS}\) themselves need to be merged which will require a renumbering. During the merge, the \(\overline\beta\) intermediates are kept sorted in hash order so that identical \(\overline\beta\) trees are deduplicated easily. The result is the a remapping of the top level list to \(\overline B_\textit{LHS} \cup \overline B_\textit{RHS}\).
    1. Finally, the checksum of the LHS and RHS is computed together and attached to the new root node.
1. A unary operation that is marked *transcendental* is immediately promoted to an intermediate \(\overline\beta\).

#### Output

The output is, for each user-defined \(\overline\gamma\), an expression tree containing a primary tree and a list of subexpressions \(\overline B_\gamma\).

### Expression tree evaluation

#### Input

A set \(\overline B_\overline\gamma\) of CSE-eliminated trees from the previous construction phase (each associated with a \(\overline\gamma\)), along with the original \(\overline \gamma\) as a function of both \(\overline\alpha \in \overline A\) and \(\overline\beta \in \overline B\).

#### Algorithm

1. At this point, we've identified subexpressions from the computation of each element of \(\overline\Gamma\) as sets \(\overline B\) that need to be computed (mapped to \(B\)). First, we need to perform a final merge as subexpressions identified in distinct \(\overline\gamma\) may overlap.
2. An issue immediately presents itself as it is possible that a subexpression of some \(\overline\gamma_k\) is located in the subexpression of a different \(\overline\gamma_l\). Thus, we perform a final N-way merge step (as when we merge Merkle trees in a binary operation) to collapse expressions one more time.
3. The map \(\overline B\rightarrow B\) is evaluated and stored to an array \(B\).
4. \(B\) is concatenated to \(A\) to form the set of runtime-values \(A\cup B\).
5. The final map \(\overline \Gamma \rightarrow \Gamma\) is evaluated and returned.

#### Output

The algorithm as outlined above has now produced the outputs as specified by \(\lambda\).

[^1]: The *indeterminate form* of a multivariate quantity substitutes variables for each known parameter to get a compile-time polynomial expression instead of a runtime valued quantity.
[^2]: ... in spite of the fact that such a lambda *is* a constant expression, but the limitations of the current state of affairs with C++ `constexpr` are plentiful.
[^3]: This is predicated on the fact that input entities are entirely disjoint (don't share any common values) at the start of the computation. At this time, this is a restriction imposed by GAL (restriction only in an optimization sense). Another way of phrasing this constraint is that all floating-point values input into the expression evaluator are assumed distinct, no matter the circumstance.
[^4]: Merkle tree allow quick \(O(1)\) comparisons of subtrees by accumulating hashes as you ascend up the tree ([wiki](https://en.wikipedia.org/wiki/Merkle_tree)). For this application, we use checksums (a CRC function) and not a hash function because we are interested primarily in fast equality checking. That is to say we aren't interested in any guarantees in the entropy of the output.
[^5]: A bloom filter is a constant-memory probabalistic data structure that supports constant time existense plausibility tests ([wiki](https://en.wikipedia.org/wiki/Bloom_filter)). If used in this application, the bloom filter implementation will be adjusted to tune for low membership and evaluation speed, not collision probability.
[^6]: Note that we do not need to repeat this procedure with the LHS and RHS arguments flipped. Existence of a common subexpression is a mutual co-property.
[^7]: Normally, using an optimization such as a bloom filter in this way would be completely overkill. It is on the table due to the poor performance of compile-time recursion.

*[CSE]: Common Subexpression Elimination
*[LHS]: Left-hand Side
*[RHS]: Reft-hand Side