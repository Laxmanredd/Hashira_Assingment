import java.io.FileReader;
import java.math.BigInteger;
import java.util.*;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

public class SecretReconstructor {

    // Large prime for finite field (2^257 - 1, suitable for 256-bit numbers)
    private static final BigInteger MOD = new BigInteger("231584178474632390847141970017375815706539969331281128078915168015826259279871");

    public static void main(String[] args) {
        try {
            // Process both test cases
            String[] testCaseFiles = {"testcase1.json", "testcase2.json"};
            List<BigInteger> secrets = new ArrayList<>();
            List<Set<Integer>> validKeysList = new ArrayList<>();

            for (String file : testCaseFiles) {
                JSONObject input = (JSONObject) new JSONParser().parse(new FileReader(file));
                JSONObject keys = (JSONObject) input.get("keys");
                long nLong = ((Number) keys.get("n")).longValue();
                long kLong = ((Number) keys.get("k")).longValue();

                // Validate n and k
                if (nLong <= 0 || kLong <= 0 || kLong > nLong) {
                    System.out.println("❌ Invalid input for " + file + ": n and k must be positive, and k must not exceed n.");
                    return;
                }
                int n = (int) nLong;
                int k = (int) kLong;

                // Parse roots
                Map<Integer, BigInteger> roots = new HashMap<>();
                for (Object key : input.keySet()) {
                    if (key.equals("keys")) continue;
                    try {
                        int x = Integer.parseInt((String) key);
                        JSONObject rootData = (JSONObject) input.get(key);
                        String baseStr = (String) rootData.get("base");
                        String value = (String) rootData.get("value");
                        int base = Integer.parseInt(baseStr);
                        BigInteger y = new BigInteger(value, base); // Decode y from given base
                        roots.put(x, y.mod(MOD)); // Apply modulo for finite field
                    } catch (Exception e) {
                        // Skip invalid roots
                    }
                }

                if (roots.size() < k) {
                    System.out.println("❌ Insufficient valid roots for " + file + ".");
                    return;
                }

                // Generate all k-combinations
                List<int[]> combinations = generateCombinations(new ArrayList<>(roots.keySet()), k);
                Map<BigInteger, List<List<Integer>>> secretMap = new HashMap<>();

                for (int[] comb : combinations) {
                    List<Integer> keys = new ArrayList<>();
                    List<BigInteger> xVals = new ArrayList<>();
                    List<BigInteger> yVals = new ArrayList<>();

                    for (int x : comb) {
                        keys.add(x);
                        xVals.add(BigInteger.valueOf(x));
                        yVals.add(roots.get(x));
                    }

                    try {
                        BigInteger[] coeffs = solvePolynomial(xVals, yVals, k);
                        BigInteger secret = coeffs[0].mod(MOD); // Constant term
                        secretMap.putIfAbsent(secret, new ArrayList<>());
                        secretMap.get(secret).add(keys);
                    } catch (ArithmeticException e) {
                        // Skip invalid combinations
                    }
                }

                if (secretMap.isEmpty()) {
                    System.out.println("❌ No valid secret could be reconstructed for " + file + ".");
                    return;
                }

                // Find most common secret
                BigInteger commonSecret = BigInteger.ZERO;
                int maxCount = 0;
                for (Map.Entry<BigInteger, List<List<Integer>>> entry : secretMap.entrySet()) {
                    if (entry.getValue().size() > maxCount) {
                        maxCount = entry.getValue().size();
                        commonSecret = entry.getKey();
                    }
                }

                // Collect valid keys
                Set<Integer> validKeys = new TreeSet<>();
                for (List<Integer> keys : secretMap.get(commonSecret)) {
                    validKeys.addAll(keys);
                }

                secrets.add(commonSecret);
                validKeysList.add(validKeys);
            }

            // Output secrets for both test cases
            System.out.println("Secrets for both test cases:");
            System.out.println("Test Case 1 Secret: " + secrets.get(0));
            System.out.println("Test Case 1 Valid Keys: " + validKeysList.get(0));
            System.out.println("Test Case 2 Secret: " + secrets.get(1));
            System.out.println("Test Case 2 Valid Keys: " + validKeysList.get(1));

        } catch (Exception e) {
            System.out.println("❌ Error processing input: " + e.getMessage());
        }
    }

    // Generate all k-combinations of items
    public static List<int[]> generateCombinations(List<Integer> items, int k) {
        List<int[]> result = new ArrayList<>();
        combine(items, 0, k, new int[k], result);
        return result;
    }

    private static void combine(List<Integer> arr, int index, int k, int[] temp, List<int[]> result) {
        if (k == 0) {
            result.add(temp.clone());
            return;
        }
        for (int i = index; i <= arr.size() - k; i++) {
            temp[temp.length - k] = arr.get(i);
            combine(arr, i + 1, k - 1, temp, result);
        }
    }

    // Solve polynomial using Gaussian elimination over finite field
    public static BigInteger[] solvePolynomial(List<BigInteger> xVals, List<BigInteger> yVals, int k) {
        BigInteger[][] A = new BigInteger[k][k];
        BigInteger[] B = new BigInteger[k];

        // Construct Vandermonde matrix
        for (int i = 0; i < k; i++) {
            BigInteger x = xVals.get(i);
            B[i] = yVals.get(i);
            for (int j = 0; j < k; j++) {
                A[i][j] = x.modPow(BigInteger.valueOf(j), MOD);
            }
        }

        // Solve Ax = B
        return gaussianElimination(A, B);
    }

    // Gaussian elimination over finite field
    public static BigInteger[] gaussianElimination(BigInteger[][] A, BigInteger[] B) {
        int n = B.length;
        BigInteger[][] augmented = new BigInteger[n][n + 1];

        // Create augmented matrix [A|B]
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmented[i][j] = A[i][j];
            }
            augmented[i][n] = B[i];
        }

        // Forward elimination
        for (int p = 0; p < n; p++) {
            // Find pivot
            int max = p;
            for (int i = p + 1; i < n; i++) {
                if (augmented[i][p].abs().compareTo(augmented[max][p].abs()) > 0) {
                    max = i;
                }
            }

            // Swap rows
            BigInteger[] temp = augmented[p];
            augmented[p] = augmented[max];
            augmented[max] = temp;

            // Check for zero pivot
            if (augmented[p][p].equals(BigInteger.ZERO)) {
                throw new ArithmeticException("Singular matrix");
            }

            // Eliminate column
            for (int i = p + 1; i < n; i++) {
                BigInteger alpha = augmented[i][p].multiply(modInverse(augmented[p][p], MOD)).mod(MOD);
                for (int j = p; j <= n; j++) {
                    augmented[i][j] = augmented[i][j].subtract(alpha.multiply(augmented[p][j])).mod(MOD);
                }
            }
        }

        // Back substitution
        BigInteger[] x = new BigInteger[n];
        for (int i = n - 1; i >= 0; i--) {
            BigInteger sum = BigInteger.ZERO;
            for (int j = i + 1; j < n; j++) {
                sum = sum.add(augmented[i][j].multiply(x[j])).mod(MOD);
            }
            x[i] = augmented[i][n].subtract(sum).multiply(modInverse(augmented[i][i], MOD)).mod(MOD);
        }

        return x;
    }

    // Modular inverse using extended Euclidean algorithm
    public static BigInteger modInverse(BigInteger a, BigInteger m) {
        BigInteger m0 = m;
        BigInteger y = BigInteger.ZERO, x = BigInteger.ONE;

        while (!a.equals(BigInteger.ZERO)) {
            BigInteger q = m.divide(a);
            BigInteger tmp = a;
            a = m.mod(a);
            m = tmp;
            BigInteger tmpX = y;
            y = x.subtract(q.multiply(y));
            x = tmpX;
        }

        if (m.compareTo(BigInteger.ONE) > 0) {
            throw new ArithmeticException("Modular inverse does not exist");
        }

        return x.mod(m0);
    }
}  