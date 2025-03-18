import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;
import java.util.Scanner;

public class GaussSeidelSolver {
    public static void main(String[] args) {
        Scanner consoleScanner = new Scanner(System.in);
        consoleScanner.useLocale(Locale.US);

        System.out.println("Выберите способ ввода данных:");
        System.out.println("1 - С клавиатуры");
        System.out.println("2 - Из файла");
        System.out.println("3 - Сгенерировать случайную матрицу");

        int choice = consoleScanner.nextInt();

        double[][] A = null;
        double[] b = null;
        int n = 0;
        double epsilon = 0.0;

        if (choice == 1) {
            System.out.println("Введите размерность системы (n <= 20):");

            n = consoleScanner.nextInt();

            A = new double[n][n];
            b = new double[n];

            System.out.println("Введите коэффициенты матрицы A построчно:");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    String token = consoleScanner.next();
                    A[i][j] = pars_double(token);
                }
            }
            System.out.println("Введите вектор свободных членов b:");
            for (int i = 0; i < n; i++) {
                String token = consoleScanner.next();
                b[i] = pars_double(token);
            }

            System.out.println("Введите требуемую точность epsilon:");
            String token = consoleScanner.next();
            epsilon = pars_double(token);

        } else if (choice == 2) {
            System.out.println("Введите путь к файлу с данными:");
            String filePath = consoleScanner.next();

            try (Scanner fileScanner = new Scanner(new File(filePath))) {
                fileScanner.useLocale(Locale.US);
                n = fileScanner.nextInt();
                A = new double[n][n];
                b = new double[n];

                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        String token = consoleScanner.next();
                        A[i][j] = pars_double(token);

                    }
                }
                for (int i = 0; i < n; i++) {
                    String token = consoleScanner.next();
                    b[i] = pars_double(token);
                }
                String token = consoleScanner.next();
                epsilon = pars_double(token);

            } catch (FileNotFoundException e) {
                System.out.println("Файл не найден: " + e.getMessage());
                return;
            }
        } else if (choice == 3) {
            System.out.println("Введите размерность системы (n <= 20):");
            n = consoleScanner.nextInt();
            A = new double[n][n];
            b = new double[n];

            Random rand = new Random();

            // Генерю матрицу A таким образом, чтобы она была диагонально преобладающей
            for (int i = 0; i < n; i++) {
                double sumOffDiag = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        A[i][j] = rand.nextDouble() * 10;
                        sumOffDiag += Math.abs(A[i][j]);
                    }
                }
                A[i][i] = sumOffDiag + rand.nextDouble() * 10 + 1;
            }

            for (int i = 0; i < n; i++) {
                b[i] = rand.nextDouble() * 100;
            }

            System.out.println("Сгенерированная матрица A:");
            for (int i = 0; i < n; i++) {
                System.out.println(Arrays.toString(A[i]));
            }
            System.out.println("Сгенерированный вектор b:");
            System.out.println(Arrays.toString(b));

            System.out.println("Введите требуемую точность epsilon:");
            epsilon = consoleScanner.nextDouble();
        } else {
            System.out.println("Неверный выбор. Программа завершена.");
            return;
        }

        boolean canReorder = reorderForDiagonalDominance(A, b);
        if (!canReorder) {
            System.out.println("Невозможно достичь диагонального преобладания для данной матрицы.");
            return;
        }

        double matrixNorm = calculateMatrixNorm(A);
        System.out.printf("Норма матрицы (максимальная сумма по строкам): %.6f%n", matrixNorm);

        GaussSeidelResult result = gaussSeidel(A, b, epsilon, 1000);

        if (!result.isConverged) {
            System.out.println("Метод Гаусса–Зейделя не сошёлся за заданное число итераций.");
        } else {
            System.out.println("Метод сошёлся за " + result.iterations + " итераций.");
            System.out.println("Решение x:");
            System.out.println(Arrays.toString(result.x));
            System.out.println("Вектор погрешностей (|x_i^(k) - x_i^(k-1)|):");
            System.out.println(Arrays.toString(result.errors));
        }
    }

    /**
     * Метод Гаусса–Зейделя.
     * Матрица системы
     * b Вектор правых частей
     * epsilon Заданная точность
     * maxIterations Максимальное число итераций
     * Результат выполнения (решение, вектор ошибок, число итераций, признак сходимости)
     */
    public static GaussSeidelResult gaussSeidel(double[][] A, double[] b, double epsilon, int maxIterations) {
        int n = A.length;
        double[] x = new double[n];      // начальное приближение (по умолчанию нули)
        double[] oldX = new double[n];   // предыдущий шаг, для оценки погрешности
        double[] errors = new double[n]; // вектор |x_i^(k) - x_i^(k-1)|

        boolean converged = false;
        int iteration = 0;

        while (iteration < maxIterations) {
            System.arraycopy(x, 0, oldX, 0, n);

            for (int i = 0; i < n; i++) {
                double sigma = 0.0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sigma += A[i][j] * x[j];
                    }
                }
                x[i] = (b[i] - sigma) / A[i][i];
            }

            double maxDiff = 0.0;
            for (int i = 0; i < n; i++) {
                errors[i] = Math.abs(x[i] - oldX[i]);
                if (errors[i] > maxDiff) {
                    maxDiff = errors[i];
                }
            }

            iteration++;
            // Проверяем критерий окончания
            if (maxDiff < epsilon) {
                converged = true;
                break;
            }
        }
        return new GaussSeidelResult(x, errors, iteration, converged);
    }

    /**
     * Пробую переставлять строки матрицы A (и соответствующие элементы b),
     * чтобы выполнить условие диагонального преобладания.
     * Возвращает true, если удалось достичь диагонального преобладания.
     */
    public static boolean reorderForDiagonalDominance(double[][] A, double[] b) {
        int n = A.length;

        for (int i = 0; i < n; i++) {
            int pivotRow = i;
            double maxValue = Math.abs(A[i][i]);
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(A[k][i]) > maxValue) {
                    pivotRow = k;
                    maxValue = Math.abs(A[k][i]);
                }
            }
            if (pivotRow != i) {
                double[] tempRow = A[i];
                A[i] = A[pivotRow];
                A[pivotRow] = tempRow;

                double tempB = b[i];
                b[i] = b[pivotRow];
                b[pivotRow] = tempB;
            }
        }

        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += Math.abs(A[i][j]);
                }
            }
            if (Math.abs(A[i][i]) <= sum) {
                return false;
            }
        }
        return true;
    }

    /**
     * Вычисление нормы матрицы (максимальная по строкам сумма модулей).
     */
    public static double calculateMatrixNorm(double[][] A) {
        double maxRowSum = 0.0;
        for (int i = 0; i < A.length; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < A[i].length; j++) {
                rowSum += Math.abs(A[i][j]);
            }
            if (rowSum > maxRowSum) {
                maxRowSum = rowSum;
            }
        }
        return maxRowSum;
    }
    public static double pars_double(String input) {
        input = input.trim();

        if (input.contains(",")) {
            input = input.replace(',', '.');
        }
        return Double.parseDouble(input);
    }

    static class GaussSeidelResult {
        double[] x;
        double[] errors;
        int iterations;
        boolean isConverged;

        public GaussSeidelResult(double[] x, double[] errors, int iterations, boolean isConverged) {
            this.x = x;
            this.errors = errors;
            this.iterations = iterations;
            this.isConverged = isConverged;
        }
    }
}
