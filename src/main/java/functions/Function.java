package functions;

public interface Function {
    double getValue(double t);

    default Function mult(final Function function) {
        return t -> Function.this.getValue(t) * function.getValue(t);
    }

    default Function add(final Function function) {
        return t -> Function.this.getValue(t) + function.getValue(t);
    }

    default Function mult(final double number) {
        return t -> Function.this.getValue(t) * number;
    }
}
