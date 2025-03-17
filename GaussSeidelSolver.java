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
            // Считываем с клавиатуры
            System.out.println("Введите размерность системы (n <= 20):");
            n = consoleScanner.nextInt();

            A = new double[n][n];
            b = new double[n];

            System.out.println("Введите коэффициенты матрицы A построчно:");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    A[i][j] = consoleScanner.nextDouble();
                }
            }
            System.out.println("Введите вектор свободных членов b:");
            for (int i = 0; i < n; i++) {
                b[i] = consoleScanner.nextDouble();
            }

            System.out.println("Введите требуемую точность epsilon:");
            epsilon = consoleScanner.nextDouble();

        } else if (choice == 2) {
            // Считываем из файла
            System.out.println("Введите путь к файлу с данными:");
            String filePath = consoleScanner.next();

            try (Scanner fileScanner = new Scanner(new File(filePath))) {
                fileScanner.useLocale(Locale.US);
                n = fileScanner.nextInt();
                A = new double[n][n];
                b = new double[n];

                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        A[i][j] = fileScanner.nextDouble();
                    }
                }
                for (int i = 0; i < n; i++) {
                    b[i] = fileScanner.nextDouble();
                }
                epsilon = fileScanner.nextDouble();

            } catch (FileNotFoundException e) {
                System.out.println("Файл не найден: " + e.getMessage());
                return;
            }
        } else if (choice == 3) {
            // Генерация случайной матрицы
            System.out.println("Введите размерность системы (n <= 20):");
            n = consoleScanner.nextInt();
            A = new double[n][n];
            b = new double[n];

            Random rand = new Random();

            // Генерируем матрицу A таким образом, чтобы она была диагонально преобладающей
            for (int i = 0; i < n; i++) {
                double sumOffDiag = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        // Генерируем случайное число от 0 до 10 для вне диагональных элементов
                        A[i][j] = rand.nextDouble() * 10;
                        sumOffDiag += Math.abs(A[i][j]);
                    }
                }
                // Диагональный элемент выбираем больше суммы модулей остальных элементов строки
                A[i][i] = sumOffDiag + rand.nextDouble() * 10 + 1;
            }

            // Генерируем случайный вектор свободных членов b (например, значения от 0 до 100)
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

        // 1) Попытка достичь диагонального преобладания
        boolean canReorder = reorderForDiagonalDominance(A, b);
        if (!canReorder) {
            System.out.println("Невозможно достичь диагонального преобладания для данной матрицы.");
            return;
        }

        // 2) Вычисляем (например) норму матрицы A (максимальная построчная сумма модулей)
        double matrixNorm = calculateMatrixNorm(A);
        System.out.printf("Норма матрицы (максимальная сумма по строкам): %.6f%n", matrixNorm);

        // 3) Запускаем метод Гаусса–Зейделя
        GaussSeidelResult result = gaussSeidel(A, b, epsilon, 1000 /*макс число итераций*/);

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
     * @param A Матрица системы
     * @param b Вектор правых частей
     * @param epsilon Заданная точность
     * @param maxIterations Максимальное число итераций
     * @return Результат выполнения (решение, вектор ошибок, число итераций, признак сходимости)
     */
    public static GaussSeidelResult gaussSeidel(double[][] A, double[] b, double epsilon, int maxIterations) {
        int n = A.length;
        double[] x = new double[n];      // начальное приближение (по умолчанию нули)
        double[] oldX = new double[n];   // предыдущий шаг, для оценки погрешности
        double[] errors = new double[n]; // вектор |x_i^(k) - x_i^(k-1)|

        boolean converged = false;
        int iteration = 0;

        while (iteration < maxIterations) {
            // Копируем x в oldX
            System.arraycopy(x, 0, oldX, 0, n);

            // Один проход по формуле Гаусса–Зейделя
            for (int i = 0; i < n; i++) {
                double sigma = 0.0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sigma += A[i][j] * x[j];
                    }
                }
                x[i] = (b[i] - sigma) / A[i][i];
            }

            // Считаем текущие ошибки
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
     * Пытаеемся переставлять строки матрицы A (и соответствующие элементы b),
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

    /**
     * Вспомогательный класс для возврата результата из метода Гаусса–Зейделя.
     */
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
