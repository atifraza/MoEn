import java.util.Comparator;

public class PairComparator implements Comparator<Pair> {
    
    @Override
    public int compare(Pair o1, Pair o2) {
        if (o1.dist < o2.dist) {
            return -1;
        } else if (o1.dist > o2.dist) {
            return 1;
        } else {
            return 0;
        }
    }
}
