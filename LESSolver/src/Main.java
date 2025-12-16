import java.util.Scanner;

public class Main{
    public static void main(){
        Scanner sc = new Scanner(System.in);

        System.out.print("Enter matrix: ");
        String input = sc.nextLine();

        double[][] A = Lessolver.parseMatrix(input);
        String result = Lessolver.solve(A);

        System.out.println(result);
    }
}
