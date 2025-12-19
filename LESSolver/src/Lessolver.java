public class Lessolver {
    static final private double precision = 1e-12;

    static String solve(double[][] Ab){
        try{
            String solutionLogs = "";

            //initial matrix log
            solutionLogs += "Initial matrix Ab:\n" + matrixToString(Ab) + "\n";
            solutionLogs += "Solving:\n\n";
            //

            int m = Ab.length, n = Ab[0].length;

            ////Gauss algorithm
            ///zeros below the main diagonal

            //logs
            solutionLogs += "0 below the main diagonal:\n";
            //

            for (int i = 0; i < m; i++) {
                //no element on the main diagonal can be equal to 0
                if (Ab[i][i] == 0){
                    for (int j = i + 1; j < m; j++){
                        if (Ab[j][i] != 0) {
                            swapRows(Ab, i, j);

                            //logs
                            solutionLogs += "Rows swaped: " + (i + 1) + " and " + (j + 1) + ":\n" + matrixToString(Ab);
                            break;
                        }
                    }
                }

                if (Ab[i][i] == 0) continue;

                for (int j = i + 1; j < m; j++){
                    //double k = - Ab[j][i] / Ab[i][i];

                    double k1 = Ab[i][i];
                    double k2 = Ab[j][i];

                    for (int w = i; w < n; w++){
                        Ab[j][w] = Ab[j][w] * k1 - Ab[i][w] * k2; //cleanPrecisionError(Ab[j][w] + k * Ab[i][w]);
                    }

                    //logs
                    //solutionLogs += "Row" + (j + 1) + " <- Row" + (j + 1) + " + " + k + " * Row " + (i + 1) + ":\n" + matrixToString(Ab);
                    solutionLogs += "Row" + (j + 1) + " <- Row" + (j + 1) + " * " + k1 + " + " + "Row" + (i + 1) + " * " + k2 + ":\n" + matrixToString(Ab);
                }
            }

            ///zeros above the main diagonal

            //logs
            solutionLogs += "\n0 above the main diagonal:\n";
            //

            for (int i = m - 1; i > 0; i--) {
                if (Ab[i][i] == 0) continue;

                for (int j = i - 1; j >= 0; j--){
                    //double k = - Ab[j][i] / Ab[i][i];

                    double k1 = Ab[i][i];
                    double k2 = Ab[j][i];
                    for (int w = n - 1; w >= 0; w--){
                        Ab[j][w] = Ab[j][w] * k1 - Ab[i][w] * k2;//cleanPrecisionError(Ab[j][w] + k * Ab[i][w]);
                    }

                    //logs
                    //solutionLogs += "Row" + (j + 1) + " <- Row" + (j + 1) + " + " + k + " * Row " + (i + 1) + ":\n" + matrixToString(Ab);
                    solutionLogs += "Row" + (j + 1) + " <- Row" + (j + 1) + " * " + k1 + " + " + "Row" + (i + 1) + " * " + k2 + ":\n" + matrixToString(Ab);
                }
            }

            ///the main diagonal equals to 1

            //logs
            solutionLogs += "\nThe main diagonal equals to 1:\n";
            //

            for (int i = 0; i < m; i++) {
                if (Ab[i][i] == 0) continue;

                double l = Ab[i][i];
                for (int w = i; w < n; w++){
                    Ab[i][w] = Ab[i][w] / l; //cleanPrecisionError(Ab[i][w] / l);
                }

                //logs
                solutionLogs += "Row" + (i + 1) + " <- Row" + (i + 1) + " * (1 / " + l + "):\n" + matrixToString(Ab);
            }

            ///Analyze the solution
            char solutionType = '1'; // '0' - no solution, '1' - there are any solutions (or exactly one)
            int rgA = 0; // rang of the matrix A
            for (int i = 0; i < m; i++) {
                if (Ab[i][i] == 0 && Ab[i][n - 1] != 0) {
                    solutionType = '0';
                    break;
                }

                if (Ab[i][i] == 0 && rgA == 0) rgA = i;
            }
            if (rgA == 0) rgA = m;  //if no 0 elements on the main diagonal were found - rang of A is equal to m (number of rows)

            if (solutionType == '0') {
                return solutionLogs + "\n------------------------------------------\nNo solution!";
            }
            else {
                String solution = "";

                //one static solution x0, number of variables (size of vector x0) = n - 1 = numberOfColsIn_Ab - numberOfColsIn_b
                double[][] x0 = new double[1][n - 1];
                for (int i = 0; i < m; i++){
                    x0[0][i] = Ab[i][n - 1];
                }

                //array of dim(ker(f)) = s = (n - 1) - rgA solutions of the A * x = 0 homogenous system of linear equations (for parameter representation)
                //n is the number of columns in matrix Ab, n - 1 is the number of columns in matrix A
                int s = (n - 1) - rgA;

                double[][] homSolutions = new double[s][];
                for (int i = rgA; i < n - 1; i++){
                    homSolutions[i - rgA] = new double[n - 1];
                    for (int j = 0; j < rgA; j++){
                        homSolutions[i - rgA][j] = - Ab[j][i];
                    }

                    homSolutions[i - rgA][i] = 1; //1 in one positon of a canonical basis vector
                }

                String variablesCap = "";
                for (int i = 1; i <= n - 1; i++){
                    variablesCap += String.format("%8s", "x" + i);
                }

                solution = "x = " + matrixToString(x0);
                for (int i = 0; i < homSolutions.length; i++){
                    solution += "+ t" + (i + 1) + " * " + matrixToString(new double[][] { homSolutions[i] });
                }

                return solutionLogs + "\n------------------------------------------\nSolution:\n" + variablesCap + "\n" + solution;
            }
        }
        catch (Exception ex){
            return "Sorry, unexpected error occurred! Try to reach the developer...\nDetails:\n" + ex.toString();
        }
    }

    static double[][] parseMatrix(String inputText){
        try {
            int numberOfOpeningBrackets = (int) inputText.chars().filter(x -> x == '[').count();
            int numberOfClosingBrackets = (int) inputText.chars().filter(x -> x == ']').count();

            if (numberOfOpeningBrackets != numberOfClosingBrackets) System.out.print("Error: Brackets mismatch!");

            //remove all spaces, tabs etc...
            inputText = inputText.replaceAll("\\s+", "");

            //number of rows
            int rows = numberOfOpeningBrackets;

            //count the number of columns
            int cols = 0;
            for (int s = 0; s < inputText.length(); s++) {
                if (Character.isDigit(inputText.charAt(s)) && (inputText.charAt(s + 1) == ']' || inputText.charAt(s + 1) == ',')) cols++;
                if (inputText.charAt(s) == ']') break;
            }

            //create a matrix
            double[][] matrix = new double[rows][cols];

            //fill matrix with values
            inputText = inputText.replaceAll("[\\[\\],]", " ").trim(); //remove brackets and commas
            String[] valuesArray = inputText.split("\\s+");

            int i = 0, j = 0;

            for (int v = 0; v < valuesArray.length; v++) {
                matrix[i][j] = Double.parseDouble(valuesArray[v]);

                if (j == cols - 1){
                    j = 0;
                    i++;
                }
                else {
                    j++;
                }
            }

            return matrix;
        }
        catch (Exception ex){
            System.out.print("Unexpected error occurred by parsing.\nDetails:\n" + ex.toString());
            return null;
        }
    }

    static void print(double[][] A) {
        System.out.print(matrixToString(A));
    }

    static String matrixToString(double[][] A){
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < A.length; i++) {
            sb.append(
                    i == 0 ? "⎡ " :
                            i == A.length - 1 ? "⎣ " : "⎢ "
            );

            for (double v : A[i]) {
                sb.append(String.format("%6.2f ", v));
            }

            sb.append(
                    i == 0 ? "⎤" :
                            i == A.length - 1 ? "⎦" : "⎥"
            );

            sb.append(System.lineSeparator());
        }

        return sb.toString();
    }

    static void swapRows(double[][] A, int row1, int row2){
        double[] tempRow = A[row1];
        A[row1] = A[row2];
        A[row2] = tempRow;
    }

//    static private double cleanPrecisionError(double value) { //clean the number when it is too small, what most likely means that it is the consequence of the finite precision of the type double. "Cleaning" means round to 0
//        return Math.abs(value) < precision ? 0.0 : value;
//    }

}
