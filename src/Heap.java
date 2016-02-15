import java.util.Arrays;

class Heap {
    private Pair[] heapArray;
    private int maxSize;
    public int currentSize;
    
    public Heap(int maxHeapSize) {
        maxSize = maxHeapSize;
        currentSize = 0;
        heapArray = new Pair[maxSize];
    }
    
    public boolean insert(Pair value) {
        if (currentSize == maxSize)
            return false;
            
        heapArray[currentSize] = value;
        cascadeUp(currentSize++);
        return true;
    }
    
    public boolean append(Pair value) {
        if (currentSize == maxSize)
            return false;
            
        heapArray[currentSize++] = value;
        return true;
    }
    
    public Pair remove() // Remove maximum value Pair
    {
        Pair root = heapArray[0];
        heapArray[0] = heapArray[--currentSize];
        cascadeDown(0);
        return root;
    }
    
    public boolean isEmpty() {
        return currentSize == 0;
    }
    
    public Pair first() {
        for (int i = 0; i < currentSize; i++)
            if (heapArray[i].dist != 0)
                return heapArray[i];
        return heapArray[0];
    }
    
    public Pair last() {
        return heapArray[currentSize - 1];
    }
    
    public Pair valueAt(int i) {
        return heapArray[i];
    }
    
    public void convertToSorted() {
        PairComparator myComparer = new PairComparator();
        Arrays.sort(heapArray, 0, currentSize, myComparer);
    }
    
    public void cascadeUp(int index) {
        int parent = (index - 1) / 2;
        Pair bottom = heapArray[index];
        while (index > 0 && heapArray[parent].dist < bottom.dist) {
            heapArray[index] = heapArray[parent];
            index = parent;
            parent = (parent - 1) / 2;
        }
        heapArray[index] = bottom;
    }
    
    public void cascadeDown(int index) {
        int largerChild;
        Pair top = heapArray[index];
        while (index < currentSize / 2) {
            int leftChild = 2 * index + 1;
            int rightChild = leftChild + 1;
            if (rightChild < currentSize
                && heapArray[leftChild].dist < heapArray[rightChild].dist)
                largerChild = rightChild;
            else
                largerChild = leftChild;
            if (top.dist >= heapArray[largerChild].dist)
                break;
            heapArray[index] = heapArray[largerChild];
            index = largerChild;
        }
        heapArray[index] = top;
    }
    
    public boolean heapIncreaseDecreaseKey(int index, Pair newValue) {
        if (index < 0 || index >= currentSize)
            return false;
        Pair oldValue = heapArray[index];
        heapArray[index] = newValue;
        if (oldValue.dist < newValue.dist)
            cascadeUp(index);
        else
            cascadeDown(index);
        return true;
    }
    
    public void displayHeap() {
        System.out.println();
        System.out.print("Elements of the Heap Array are : ");
        for (int m = 0; m < currentSize; m++)
            if (heapArray[m].dist != 0.0)
                System.out.print(heapArray[m].dist + " ");
            else
                System.out.print("-- ");
        System.out.println();
        int emptyLeaf = 32;
        int itemsPerRow = 1;
        int column = 0;
        int j = 0;
        String separator = "...............................";
        System.out.println(separator + separator);
        while (currentSize > 0) {
            if (column == 0)
                for (int k = 0; k < emptyLeaf; k++)
                    System.out.print(' ');
            System.out.print(heapArray[j].dist);
            
            if (++j == currentSize)
                break;
            if (++column == itemsPerRow) {
                emptyLeaf /= 2;
                itemsPerRow *= 2;
                column = 0;
                System.out.println();
            } else
                for (int k = 0; k < emptyLeaf * 2 - 2; k++)
                    System.out.print(' ');
        }
        System.out.println("\n" + separator + separator);
    }
}
