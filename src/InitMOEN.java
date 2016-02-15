import java.io.BufferedReader;
import java.io.FileReader;
import java.time.Duration;
import java.time.Instant;

public class InitMOEN {
    static int cap;
    static double ccc;
    static int VERBOSE = 0;
    static int WW = 0;
    static int K = 6;
    
    public static void main(String[] args) {
        double[] TS;
//        double[] xy;
        double[] x;
        double[] x2;
//        double[] prevXY;
        
        double[][] bsf;
        long[][] loc1;
        long[][] loc2;
        
        double[] maX;
        
        Instant Start = Instant.now();
        Duration Elapsed;
        int SIZE = Integer.parseInt(args[1]);
        
        int n = 1;
        int minLength = Integer.parseInt(args[2]);
        int maxLength = Integer.parseInt(args[3]);
        K = 50;
        if (args.length > 4)
            K = Integer.parseInt(args[4]);
        ccc = 0.1;
        if (args.length > 5)
            ccc = Double.parseDouble(args[5]);
        WW = SIZE / 2;
        cap = SIZE;// *(int)minLength;
        
        try {
            BufferedReader sr = new BufferedReader(new FileReader(args[0]));
            TS = new double[SIZE];
            x2 = new double[SIZE];
            x = new double[SIZE];
//            xy = new double[SIZE];
//            prevXY = new double[SIZE];
            
            bsf = new double[maxLength][K];
            loc1 = new long[maxLength][K];
            loc2 = new long[maxLength][K];
            
            maX = new double[maxLength];
            
            TS[0] = 0;
            x[0] = 0;
            x2[0] = 0;
            String line;
            
            while ((line = sr.readLine()) != null && n < SIZE - 2) {
                double d = Double.parseDouble(line);
                TS[n] = d;
                x[n] = d + x[n - 1];
                x2[n] = d * d + x2[n - 1];
                
                n++;
            }
            sr.close();
            
            if (VERBOSE == 1)
                System.out.println("The size of the data is : " + n);
                
            // Initialize bsf
            for (int j = 0; j < maxLength; j++)
                for (int k = 0; k < K; k++)
                    bsf[j][k] = Double.MAX_VALUE;
                    
            int count = 0;
            int mistakeCount = 0;
            
            /* The new algorithm for the length independent motif */
            
            Start = Instant.now();
            count = 0;
            
            //// Start: maX computation block
            
            for (int j = 0; j < maxLength; j++) {
                maX[j] = 0;
            }
            
            for (int i = 1; i < n - maxLength; i++) {
                int mLoc = -1;
                int mnLoc = -1;
                double runMax = 0;
                double runMin = 999999999999999.099;
                
                for (int j = minLength; j < maxLength; j++) {
                    
                    double sumY = x[i + j - 1] - x[i - 1];
                    double meanY = sumY / j;
                    double sumY2 = x2[i + j - 1] - x2[i - 1];
                    double sigmaY = Math.sqrt((sumY2 / j) - (meanY * meanY));
                    
                    if (j == minLength) {
                        runMax = 0;
                        runMin = 999999999999999.099;
                        for (int k = 1; k <= minLength; k++) {
                            double X = (TS[i + k - 1] - meanY) / sigmaY;
                            if (X > runMax) {
                                runMax = X;
                                mLoc = i + k - 1;
                            }
                            if (X < runMin) {
                                runMin = X;
                                mnLoc = i + k - 1;
                            }
                            
                        }
                        
                    } else {
                        runMax = (TS[mLoc] - meanY) / sigmaY;
                        runMin = (TS[mnLoc] - meanY) / sigmaY;
                        double Y = (TS[i + j - 1] - meanY) / sigmaY;
                        if (runMax < Y) {
                            runMax = Y;
                            mLoc = i + j - 1;
                        }
                        if (runMin > Y) {
                            runMin = Y;
                            mnLoc = i + j - 1;
                        }
                        
                    }
                    
                    if (runMax > maX[j])
                        maX[j] = runMax;
                    // if (-runMax*runMin > maX[j])
                    // maX[j] = -runMax * runMin;
                    if (-runMin > maX[j])
                        maX[j] = -runMin;
                        
                }
                
                // System.out.println(i+ " " + maX[128]);
            }
            
            Elapsed = Duration.between(Instant.now(), Start);
            if (VERBOSE == 1)
                System.out.println("Time Elapsed till Max Computation : {0} ms"
                                   + Elapsed.getSeconds());
                                   
            /// END : maX computation block
            ///
            
            Heap List = null;
            double firstLB = -1;
            for (int j = minLength; j < maxLength; j++) {
                
                // System.out.println("LENGTH " + j);
                if (List == null || List.isEmpty() == true) {
                    // List = FindMotifNewSpace(TS, SIZE, j);
                    List = FindMotif(TS, SIZE, j);
                    List.convertToSorted();
                    
                    double dis = List.last().dist;
                    double z = maX[j];
                    double LB = j / (j - 1) + z * z * (j - 1) / (j * j);
                    firstLB = dis / Math.sqrt(LB);
                    Elapsed = Duration.between(Instant.now(), Start);
                    
                    int kk = 0;
                    Pair curr = (Pair) List.valueAt(0);
                    loc1[j][kk] = curr.loc1;
                    loc2[j][kk] = curr.loc2;
                    bsf[j][kk] = curr.dist;
                    kk++;
                    for (int nk = 1; nk < List.currentSize && kk < K; nk++) {
                        curr = (Pair) List.valueAt(nk);
                        int sign = 0;
                        for (int mk = 0; mk < kk && kk < K; mk++) {
                            Pair test = (Pair) List.valueAt(mk);
                            // long d1 = j - Math.Abs(curr.loc2 -
                            // test.loc2);
                            // long d2 = j - Math.Abs(curr.loc1 -
                            // test.loc1);
                            if (isCovering(curr, (int) j, test, (int) j))
                            // if (d1 < 0.2*j || d2 > 0.2*j)
                            {
                                sign = 1;
                                break;
                            }
                        }
                        if (sign == 0) {
                            loc1[j][kk] = curr.loc1;
                            loc2[j][kk] = curr.loc2;
                            bsf[j][kk] = curr.dist;
                            kk++;
                        }
                    }
                    
                    if (VERBOSE == 1)
                        System.out.println("Time Elapsed till First Motif Computation : {0} ms"
                                           + Elapsed.getSeconds());
                    count++;
                    // printList(List);
                } else {
                    
                    double LB = 0, best = 0;
                    double dis = List.last().dist;
                    double z = maX[j];
                    LB = j / (j - 1) + z * z * (j - 1) / (j * j);
                    LB = dis / Math.sqrt(LB);
                    
                    if (VERBOSE == 1)
                        
                        System.out.println(firstLB + " FIRSTLB " + LB
                                           + " : LB   " + dis + " : dists\n");
                    Heap nextList = new Heap(List.currentSize);
                    
                    for (int k = 0; k < List.currentSize; k++) {
                        
                        Pair curr = (Pair) List.valueAt(k);
                        if (curr.loc1 + j >= n || curr.loc2 + j >= n
                            || curr.loc1 == 0 || curr.loc2 == 0)
                            continue;
                        curr.dist = distance(TS, curr.loc1, curr.loc2, j);
                        
                        nextList.append(curr);
                        
                        // System.out.println(k);
                    }
                    
                    nextList.convertToSorted();
                    
                    int kk = 0;
                    Pair curr1 = (Pair) List.valueAt(0);
                    best = curr1.dist;
                    kk++;
                    for (int nk = 1; nk < List.currentSize && kk < K; nk++) {
                        curr1 = (Pair) List.valueAt(nk);
                        int sign = 0;
                        for (int mk = 0; mk < kk && kk < K; mk++) {
                            Pair test = (Pair) List.valueAt(mk);
                            // long d1 = j - Math.Abs(curr1.loc2 -
                            // test.loc2);
                            // long d2 = j - Math.Abs(curr1.loc1 -
                            // test.loc1);
                            if (isCovering(curr1, (int) j, test, (int) j)) {
                                sign = 1;
                                break;
                            }
                        }
                        if (sign == 0) {
                            best = curr1.dist;
                            kk++;
                        }
                    }
                    
                    if (best >= firstLB) {
                        
                        List = FindMotif(TS, SIZE, j);
                        List.convertToSorted();
                        dis = List.last().dist;
                        z = maX[j];
                        LB = j / (j - 1) + z * z * (j - 1) / (j * j);
                        firstLB = dis / Math.sqrt(LB);
                        if (VERBOSE == 1)
                            System.out.println("\n\nFIRST LB " + firstLB
                                               + "\n\n");
                                               
                        count++;
                        
                    } else if (best < firstLB) {
                        
                        // Cut down the list to maintain order. ||
                        // (double)List.GetKey(i) > LB
                        
                        List = nextList;
                        z = maX[j];
                        LB = j / (j - 1) + z * z * (j - 1) / (j * j);
                        firstLB = firstLB / Math.sqrt(LB);
                        
                    }
                    kk = 0;
                    curr1 = (Pair) List.valueAt(0);
                    loc1[j][kk] = curr1.loc1;
                    loc2[j][kk] = curr1.loc2;
                    bsf[j][kk] = curr1.dist;
                    kk++;
                    for (int nk = 1; nk < List.currentSize && kk < K; nk++) {
                        curr1 = (Pair) List.valueAt(nk);
                        int sign = 0;
                        for (int mk = 0; mk < kk && kk < K; mk++) {
                            Pair test = (Pair) List.valueAt(mk);
                            // long d1 = j - Math.Abs(curr1.loc2 -
                            // test.loc2);
                            // long d2 = j - Math.Abs(curr1.loc1 -
                            // test.loc1);
                            // if (d1 > 0.2 * j || d2 > 0.2 * j)
                            if (isCovering(curr1, (int) j, test, (int) j)) {
                                sign = 1;
                                break;
                            }
                        }
                        if (sign == 0) {
                            
                            loc1[j][kk] = curr1.loc1;
                            loc2[j][kk] = curr1.loc2;
                            bsf[j][kk] = curr1.dist;
                            kk++;
                            
                        }
                    }
                    
                    Pair P = (Pair) List.valueAt(K - 1);
                    if (VERBOSE == 1) {
                        System.out.println("Correct Numbers " + bsf[j][K - 1]
                                           + "\t" + loc1[j][K - 1] + "\t"
                                           + loc2[j][K - 1]);
                        System.out.println(P.dist + "\t" + P.loc1 + "\t"
                                           + P.loc2 + "\t" + j + "\t"
                                           + List.currentSize);
                    }
                    
                }
                
            }
            
            Elapsed = Duration.between(Instant.now(), Start);
            if (VERBOSE == 1) {
                System.out.println("Count of FindMotif calls : " + count);
                System.out.println("Count of Guesses : "
                                   + (maxLength - minLength - count));
                System.out.println("mistakeCount " + mistakeCount);
                System.out.println("Time Elapsed: {0} ms"
                                   + Elapsed.getSeconds());
            }
            
            int[][] mark = new int[maxLength][K];
            
            for (int i = minLength; i < maxLength; i++) {
                
//                double t = Math.floor(i * 0.2);
                for (int kk = 0; kk < K; kk++) {
                    mark[i][kk] = 0;
                    int j = i;
                    for (int kj = kk + 1; kj < K; kj++)
                        if (isCovering(loc1[i][kk], loc2[i][kk], (int) j,
                                       loc1[i][kj], loc2[i][kj], (int) j))
                        // if (loc1[i, kk] - loc1[j, kj] <= j - i + t
                        // && loc2[i, kk] - loc2[j, kj] <= j - i + t
                        // && loc1[i, kk] - loc1[j, kj] >= 0 - t &&
                        // loc2[i, kk] - loc2[j, kj] >= 0 - t)
                        {
                            mark[i][kk] = 1;
                            break;
                        }
                    for (j = i + 1; j < maxLength; j++) {
                        for (int kj = 0; kj < K; kj++)
                            // if (isCovering(loc1[j, kj], loc2[j,
                            // kj], (int)j,loc1[i, kk], loc2[i, kk],
                            // (int)i) && bsf[j,kj] < Math.sqrt(0.5*j)
                            // )
                            if (isCovering(loc1[j][kj], loc2[j][kj], (int) j,
                                           loc1[i][kk], loc2[i][kk], (int) i))
                            // if (loc1[i,kk] - loc1[j,kj] <= j - i +
                            // t && loc2[i,kk] - loc2[j,kj] <= j - i +
                            // t && loc1[i,kk] - loc1[j,kj] >= 0 - t
                            // && loc2[i,kk] - loc2[j,kj] >= 0 - t)
                            {
                                mark[i][kk] = 1;
                                break;
                                
                            }
                            
                    }
                }
            }
            
            for (int i = minLength; i < maxLength; i++)
                for (int kk = 0; kk < K; kk++)
                    if (mark[i][kk] == 0)
                        System.out.println(i + "," + loc1[i][kk] + ","
                                           + loc2[i][kk] + "," + bsf[i][kk]);
                                           
        } catch (Exception e) {
            System.out.println("The file could not be read:");
            System.out.println(e);
        }
        
    }
    
    // static void printList(SortedList L)
    // {
    // for (int k = 0; k < L.Count; k++)
    // {
    // Pair PP = (Pair)L.GetByIndex(k);
    // System.out.println((double)L.GetKey(k) + "\t" + PP.loc1 + "\t"
    // + PP.loc2);
    //
    // }
    //
    // }
    
    static Heap FindMotif(double[] TS, int SIZE, int LENGTH) {
        
        try {
            double[] xy;
            double[] x;
            double[] x2;
            double[] prevXY;
            
            double bsf;
            int loc1 = -1;
            int loc2 = -1;
            
            x2 = new double[SIZE];
            x = new double[SIZE];
            xy = new double[SIZE];
            prevXY = new double[SIZE];
            
            Heap List = new Heap(cap);
            
            TS[0] = 0;
            x[0] = 0;
            x2[0] = 0;
            
            int n = 1;
            
            while (n < SIZE - 2) {
                
                double d = TS[n];
                x[n] = d + x[n - 1];
                x2[n] = d * d + x2[n - 1];
                n++;
            }
            
            // System.out.println("The size of the data is : " + n);
            
            // Initialize bsf
            bsf = Double.MAX_VALUE;
            
            // compute the first dot product. Can be done quickly by
            // FFT
            xy[0] = 0;
            for (int i = 1; i <= n - LENGTH; i++) {
                xy[i] = 0;
                for (int j = 0; j < LENGTH; j++) {
                    xy[i] += TS[i + j] * TS[j + 1];
                }
                prevXY[i] = xy[i];
            }
            
            // for (long i = 1; i < n - WW - LENGTH; i++)
            for (int i = 1; i < n - LENGTH; i++) {
                // process dot product
                
                if (i > 1) {
                    for (int k = i + LENGTH; k <= n - LENGTH; k++) {
                        xy[k] = prevXY[k - 1] - TS[k - 1] * TS[i - 1]
                                + TS[k + LENGTH - 1] * TS[i + LENGTH - 1];
                        prevXY[k - 1] = xy[k - 1];
                    }
                    prevXY[n - LENGTH] = x[n - LENGTH];
                    
                }
                
                // I have a query at this point which is TS[i:i+j]
                // Find the nearest neighbor
                // for (long k = i + WW + LENGTH; k < n - LENGTH; k++)
                for (int k = i + LENGTH; k < n - LENGTH; k++) {
                    
                    int j = LENGTH;
                    // Get the necessary statistics
                    double sumX = x[k + j - 1] - x[k - 1];
                    double meanX = sumX / j;
                    double sumX2 = x2[k + j - 1] - x2[k - 1];
                    double sigmaX = Math.sqrt((sumX2 / j) - (meanX * meanX));
                    double sumY = x[i + j - 1] - x[i - 1];
                    double meanY = sumY / j;
                    double sumY2 = x2[i + j - 1] - x2[i - 1];
                    double sigmaY = Math.sqrt((sumY2 / j) - (meanY * meanY));
                    
                    // Compute Distance
                    double corr = (xy[k] - (j * meanX * meanY))
                                  / (j * sigmaX * sigmaY);
                    double dis = Math.sqrt(2 * j * (1 - (float) corr));
                    
                    // System.out.println(dis);
                    if (dis < bsf) {
                        bsf = dis;
                        loc1 = i;
                        loc2 = k;
                    }
                    
                    Pair P = new Pair();
                    P.loc1 = i;
                    P.loc2 = k;
                    P.dist = dis;
                    
                    // Insert in the sorted list
                    if (List.insert(P) == false) {
                        List.heapIncreaseDecreaseKey(0, P);
                        
                    }
                    
                }
            }
            
            System.out.println(bsf + "\t" + loc1 + "\t" + loc2 + "\t" + LENGTH);
            return List;
            
        } catch (Exception e) {
            System.out.println("ERROR in the findMotif function!!!");
            System.out.println(e);
        }
        
        return null;
    }
    
    static double distanceEAB(double[][] data, int i, int k, int length,
                              double bsf) // double bsf =
                                          // 999999999.999
    {
        double sum = 0;
        double bsf2 = bsf * bsf;
        int j;
        for (j = 0; j < length && sum <= bsf2; j++) {
            sum += (data[i][j] - data[k][j]) * (data[i][j] - data[k][j]);
            
        }
        return Math.sqrt(sum);
    }
    
    static boolean isCovering(Pair large, int largeLength, Pair small,
                              int smallLength) {
        // assumes loc1 < loc2
        double c = ccc;
        if (large.loc1 - small.loc1 >= 0
            && large.loc1 - small.loc1 < (1 - c) * smallLength)
            return true;
        else if (largeLength - small.loc1 + large.loc1 > c * smallLength
                 && largeLength - small.loc1 + large.loc1 <= largeLength)
            return true;
        else if (large.loc2 - small.loc2 >= 0
                 && large.loc2 - small.loc2 < (1 - c) * smallLength)
            return true;
        else if (largeLength - small.loc2 + large.loc2 > c * smallLength
                 && largeLength - small.loc2 + large.loc2 <= largeLength)
            return true;
        else if (large.loc1 - small.loc2 >= 0
                 && large.loc1 - small.loc2 < (1 - c) * smallLength)
            return true;
        else if (largeLength - small.loc2 + large.loc1 > c * smallLength
                 && largeLength - small.loc2 + large.loc1 <= largeLength)
            return true;
        else if (large.loc1 - small.loc2 >= 0
                 && large.loc1 - small.loc2 < (1 - c) * smallLength)
            return true;
        else if (largeLength - small.loc2 + large.loc1 > c * smallLength
                 && largeLength - small.loc2 + large.loc1 <= largeLength)
            return true;
        else
            return false;
    }
    
    static boolean isCovering(long largeLoc1, long largeLoc2, int largeLength,
                              long smallLoc1, long smallLoc2, int smallLength) {
        // assumes loc1 < loc2
        
        double c = ccc;
        if (largeLoc1 - smallLoc1 >= 0
            && largeLoc1 - smallLoc1 < (1 - c) * smallLength)
            return true;
        else if (largeLength - smallLoc1 + largeLoc1 > c * smallLength
                 && largeLength - smallLoc1 + largeLoc1 <= largeLength)
            return true;
        else if (largeLoc2 - smallLoc2 >= 0
                 && largeLoc2 - smallLoc2 < (1 - c) * smallLength)
            return true;
        else if (largeLength - smallLoc2 + largeLoc2 > c * smallLength
                 && largeLength - smallLoc2 + largeLoc2 <= largeLength)
            return true;
        else if (largeLoc1 - smallLoc2 >= 0
                 && largeLoc1 - smallLoc2 < (1 - c) * smallLength)
            return true;
        else if (largeLength - smallLoc2 + largeLoc1 > c * smallLength
                 && largeLength - smallLoc2 + largeLoc1 <= largeLength)
            return true;
        else if (largeLoc1 - smallLoc2 >= 0
                 && largeLoc1 - smallLoc1 < (1 - c) * smallLength)
            return true;
        else if (largeLength - smallLoc2 + largeLoc1 > c * smallLength
                 && largeLength - smallLoc2 + largeLoc1 <= largeLength)
            return true;
            
        else
            return false;
    }
    
    static double distance(double[] TS, int i, int k, int length) {
        double xy = 0, x = 0, y = 0, x2 = 0, y2 = 0;
        for (int ii = 0; ii < length; ii++) {
            xy += TS[ii + i] * TS[ii + k];
            x += TS[ii + i];
            x2 += TS[ii + i] * TS[ii + i];
            
            y += TS[ii + k];
            y2 += TS[ii + k] * TS[ii + k];
            
        }
        
        double meanX = x / length;
        double sigmaX = Math.sqrt((x2 / length) - (meanX * meanX));
        double meanY = y / length;
        double sigmaY = Math.sqrt((y2 / length) - (meanY * meanY));
        
        // Compute Distance
        double corr = (xy - (length * meanX * meanY))
                      / (length * sigmaX * sigmaY);
        double dis = Math.sqrt(2 * length * (1 - (float) corr));
        return dis;
        
    }
}
