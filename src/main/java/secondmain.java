import functions.Function;

public class secondmain {
    public static void main(String[] args) {
        System.out.println(((Function) t -> t + 1).add(t -> 2 * t).getValue(1));
    }
}
