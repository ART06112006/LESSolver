# Linear Equations Systems Solver (LESSolver)

**LESSolver** is a specialized tool designed to solve systems of linear equations using the **Gaussian elimination algorithm**. It provides a clear, step-by-step breakdown of the computation process and handles various types of mathematical outcomes.

## Core Functionality

* **Gaussian Elimination:** Implements the classic row reduction method to find solutions.
* **Comprehensive Analysis:** Accurately identifies:
    * A unique solution.
    * Infinitely many solutions (parameters are represented as "t"`s with appropriated indices).
    * No solution (inconsistent systems).
* **High Precision:** Operates with `double` precision types.
* **Step-by-Step Visualization:** Displays the transformation of the matrix at each stage of the calculation.

---

## Input Format

To use the solver, you must provide the **augmented matrix** $(A|b)$ in a single line. Each row of the matrix should be enclosed in square brackets `[]`, and rows must be separated by commas.

### Input Syntax:
`[row1], [row2], ..., [rowN]`

**Example:**
To solve a system represented by:

x1 + 2 * x2 - x3 = 4

3 * x2 + 2 * x3 = 5

2 * x1 - x2 + x3 = 3



You should enter:
```text
[1, 2, -1, 4], [0, 3, 2, 5], [2, -1, 1, 3]
