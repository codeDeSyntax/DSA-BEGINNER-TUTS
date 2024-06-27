import java.util.InputMismatchException;
import java.util.Scanner;

public class AlgorithmFramework {

    public static int bruteForceMaxSubarraySum(int[] arr) {
        int maxSum = Integer.MIN_VALUE;
        int n = arr.length;

        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                int sum = 0;
                for (int k = i; k <= j; k++) {
                    sum += arr[k];
                }
                maxSum = Math.max(maxSum, sum);
            }
        }

        return maxSum;
    }

    public static int bruteForceStringMatching(String text, String pattern) {
        int n = text.length();
        int m = pattern.length();

        for (int i = 0; i <= n - m; i++) {
            int j;
            for (j = 0; j < m; j++) {
                if (text.charAt(i + j) != pattern.charAt(j)) {
                    break;
                }
            }
            if (j == m) {
                return i; // Match found at index i
            }
        }

        return -1; // No match found
    }

    public static double bruteForceClosestPair(double[][] points) {
        int n = points.length;
        double minDist = Double.MAX_VALUE;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double dist = Math.sqrt(Math.pow(points[i][0] - points[j][0], 2) + Math.pow(points[i][1] - points[j][1], 2));
                minDist = Math.min(minDist, dist);
            }
        }

        return minDist;
    }

    public static int tsp(int[][] graph, boolean[] visited, int currPos, int n, int count, int cost, int ans) {
        if (count == n && graph[currPos][0] > 0) {
            return Math.min(ans, cost + graph[currPos][0]);
        }

        for (int i = 0; i < n; i++) {
            if (!visited[i] && graph[currPos][i] > 0) {
                visited[i] = true;
                ans = tsp(graph, visited, i, n, count + 1, cost + graph[currPos][i], ans);
                visited[i] = false;
            }
        }

        return ans;
    }

    public static int exhaustiveKnapsack(int[] weights, int[] values, int n, int capacity) {
        if (n == 0 || capacity == 0) {
            return 0;
        }

        if (weights[n - 1] > capacity) {
            return exhaustiveKnapsack(weights, values, n - 1, capacity);
        } else {
            return Math.max(values[n - 1] + exhaustiveKnapsack(weights, values, n - 1, capacity - weights[n - 1]),
                    exhaustiveKnapsack(weights, values, n - 1, capacity));
        }
    }

    public static void generatePermutations(String str, String ans) {
        if (str.length() == 0) {
            System.out.println(ans + " ");
            return;
        }

        for (int i = 0; i < str.length(); i++) {
            char ch = str.charAt(i);
            String ros = str.substring(0, i) + str.substring(i + 1);
            generatePermutations(ros, ans + ch);
        }
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        while (true) {
            System.out.println("Choose the type of algorithm to run:");
            System.out.println("1. Brute Force Algorithms");
            System.out.println("2. Exhaustive Search Algorithms");
            System.out.println("3. Exit");
            int algorithmType = getIntInput(scanner);

            if (algorithmType == 3) {
                break;
            }

            switch (algorithmType) {
                case 1:
                    handleBruteForceAlgorithms(scanner);
                    break;

                case 2:
                    handleExhaustiveSearchAlgorithms(scanner);
                    break;

                default:
                    System.out.println("Invalid choice. Please choose 1, 2, or 3.");
                    break;
            }
        }

        scanner.close();
        System.out.println("Exited the program.");
    }

    private static void handleBruteForceAlgorithms(Scanner scanner) {
        System.out.println("Choose the Brute Force algorithm:");
        System.out.println("1. Maximum Subarray Sum");
        System.out.println("2. String Matching");
        System.out.println("3. Closest Pair of Points");
        int bruteChoice = getIntInput(scanner);

        switch (bruteChoice) {
            case 1:
                System.out.println("Enter the size of the array:");
                int n = getIntInput(scanner);
                int[] arr = new int[n];
                System.out.println("Enter the elements of the array:");
                for (int i = 0; i < n; i++) {
                    arr[i] = getIntInput(scanner);
                }

                long startTime = System.nanoTime();
                int maxSum = bruteForceMaxSubarraySum(arr);
                long endTime = System.nanoTime();

                System.out.println("Maximum subarray sum is: " + maxSum);
                System.out.println("Running time: " + (endTime - startTime) + " nanoseconds");
                System.out.println("Time Complexity: O(n^3)");
                break;

            case 2:
                System.out.println("Enter the text:");
                String text = getStringInput(scanner);
                System.out.println("Enter the pattern:");
                String pattern = getStringInput(scanner);

                long startTime2 = System.nanoTime();
                int matchIndex = bruteForceStringMatching(text, pattern);
                long endTime2 = System.nanoTime();

                if (matchIndex != -1) {
                    System.out.println("Pattern found at index: " + matchIndex);
                } else {
                    System.out.println("Pattern not found.");
                }
                System.out.println("Running time: " + (endTime2 - startTime2) + " nanoseconds");
                System.out.println("Time Complexity: O(n*m)");
                break;

            case 3:
                System.out.println("Enter the number of points:");
                int pointsCount = getIntInput(scanner);
                double[][] points = new double[pointsCount][2];
                System.out.println("Enter the coordinates of the points (x y):");
                for (int i = 0; i < pointsCount; i++) {
                    points[i][0] = getDoubleInput(scanner);
                    points[i][1] = getDoubleInput(scanner);
                }

                long startTime3 = System.nanoTime();
                double minDist = bruteForceClosestPair(points);
                long endTime3 = System.nanoTime();

                System.out.println("Minimum distance is: " + minDist);
                System.out.println("Running time: " + (endTime3 - startTime3) + " nanoseconds");
                System.out.println("Time Complexity: O(n^2)");
                break;

            default:
                System.out.println("Invalid choice.");
                break;
        }
    }

    private static void handleExhaustiveSearchAlgorithms(Scanner scanner) {
        System.out.println("Choose the Exhaustive Search algorithm:");
        System.out.println("1. Traveling Salesman Problem (TSP)");
        System.out.println("2. Knapsack Problem");
        System.out.println("3. Permutations");
        int exhaustiveChoice = getIntInput(scanner);

        switch (exhaustiveChoice) {
            case 1:
                System.out.println("Enter the number of cities:");
                int cities = getIntInput(scanner);
                int[][] graph = new int[cities][cities];

                System.out.println("Enter the cost matrix:");
                for (int i = 0; i < cities; i++) {
                    for (int j = 0; j < cities; j++) {
                        graph[i][j] = getIntInput(scanner);
                    }
                }

                boolean[] visited = new boolean[cities];
                visited[0] = true;

                long startTime = System.nanoTime();
                int result = tsp(graph, visited, 0, cities, 1, 0, Integer.MAX_VALUE);
                long endTime = System.nanoTime();

                System.out.println("Minimum cost of traveling is: " + result);
                System.out.println("Running time: " + (endTime - startTime) + " nanoseconds");
                System.out.println("Time Complexity: O(n!)");
                break;

            case 2:
                System.out.println("Enter the number of items:");
                int items = getIntInput(scanner);
                int[] weights = new int[items];
                int[] values = new int[items];

                System.out.println("Enter the weights of the items separated by white spaces:");
                for (int i = 0; i < items; i++) {
                    weights[i] = getIntInput(scanner);
                }

                System.out.println("Enter the values of the items (white space separation):");
                for (int i = 0; i < items; i++) {
                    values[i] = getIntInput(scanner);
                }

                System.out.println("Enter the capacity of the knapsack:");
                int capacity = getIntInput(scanner);

                long startTime2 = System.nanoTime();
                int maxValue = exhaustiveKnapsack(weights, values, items, capacity);
                long endTime2 = System.nanoTime();

                System.out.println("Maximum value in knapsack is: " + maxValue);
                System.out.println("Running time: " + (endTime2 - startTime2) + " nanoseconds");
                System.out.println("Time Complexity: O(2^n)");
                break;

            case 3:
                System.out.println("Enter the string:");
                String str = getStringInput(scanner);

                long startTime3 = System.nanoTime();
                generatePermutations(str, "");
                long endTime3 = System.nanoTime();

                System.out.println("Running time: " + (endTime3 - startTime3) + " nanoseconds");
                System.out.println("Time Complexity: O(n!)");
                break;

            default:
                System.out.println("Invalid choice.");
                break;
        }
    }

    private static int getIntInput(Scanner scanner) {
        while (true) {
            try {
                return scanner.nextInt();
            } catch (InputMismatchException e) {
                System.out.println("Invalid input. Please enter an integer.");
                scanner.next(); // Clear the invalid input
            }
        }
    }

    private static double getDoubleInput(Scanner scanner) {
        while (true) {
            try {
                return scanner.nextDouble();
            } catch (InputMismatchException e) {
                System.out.println("Invalid input. Please enter a double.");
                scanner.next(); // Clear the invalid input
            }
        }
    }

    private static String getStringInput(Scanner scanner) {
        while (true) {
            try {
                return scanner.next();
            } catch (InputMismatchException e) {
                System.out.println("Invalid input. Please enter a string.");
                scanner.next(); // Clear the invalid input
            }
        }
    }
}
