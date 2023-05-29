#include <vector>
#include <cstdio>
#include <algorithm>
#include <ctime>

class Log {
public:
    static void print(const char *strategy, long long int start, long long int end) {
        printf("%s: %f seconds\n", strategy, (double) (end - start) / 1000000);
    }

    static void print(std::vector<std::vector<int> > result) {
        int n = (int) result.size();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%d ", result[i][j]);
            }
            printf("\n");
        }
    }
};

class MatrixGenerator {
private:
    int n;
    std::vector<std::vector<int> > matrix;

protected:
    int generateNumber() {
        return std::rand() % 9 + 1;
    }

public:
    MatrixGenerator(int n) {
        this->n = n;
    }

    void generate() {
        this->matrix.resize(this->n, std::vector<int>(this->n, 0));

        for (int i = 0; i < this->n; i++) {
            for (int j = 0; j < this->n; j++) {
                this->matrix[i][j] = this->generateNumber();
            }
        }
    }

    std::vector<std::vector<int> > getGeneratedMatrix() {
        return this->matrix;
    }
};

class Strategy {
public:
    virtual ~Strategy() {}

    virtual void multiply(std::vector<std::vector<int> > a, std::vector<std::vector<int> > b) {
        printf("%s\n", "Please, implement this default Strategy.");
    }
};

class IJKStrategy : public Strategy {
public:
    void multiply(std::vector<std::vector<int> > a, std::vector<std::vector<int> > b) {
        int n = (int) a.size();
        std::vector<std::vector<int> > result(n, std::vector<int>(n, 0));

        long long int start = clock();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        long long int end = clock();
        Log::print("ijk", start, end);
    }
};

class IKJStrategy : public Strategy {
public:
    void multiply(std::vector<std::vector<int> > a, std::vector<std::vector<int> > b) {
        int n = (int) a.size();
        std::vector<std::vector<int> > result(n, std::vector<int>(n, 0));

        long long int start = clock();
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                for (int j = 0; j < n; j++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        long long int end = clock();
        Log::print("ikj", start, end);
    }
};

class JIKStrategy : public Strategy {
public:
    void multiply(std::vector<std::vector<int> > a, std::vector<std::vector<int> > b) {
        int n = (int) a.size();
        std::vector<std::vector<int> > result(n, std::vector<int>(n, 0));

        long long int start = clock();
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                for (int k = 0; k < n; k++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        long long int end = clock();
        Log::print("jik", start, end);
    }
};

class JKIStrategy : public Strategy {
public:
    void multiply(std::vector<std::vector<int> > a, std::vector<std::vector<int> > b) {
        int n = (int) a.size();
        std::vector<std::vector<int> > result(n, std::vector<int>(n, 0));

        long long int start = clock();
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                for (int i = 0; i < n; i++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        long long int end = clock();
        Log::print("jki", start, end);
    }
};

class KIJStrategy : public Strategy {
public:
    void multiply(std::vector<std::vector<int> > a, std::vector<std::vector<int> > b) {
        int n = (int) a.size();
        std::vector<std::vector<int> > result(n, std::vector<int>(n, 0));
        long long int start = clock();
        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        long long int end = clock();
        Log::print("kij", start, end);
    }
};

class KJIStrategy : public Strategy {
public:
    void multiply(std::vector<std::vector<int> > a, std::vector<std::vector<int> > b) {
        int n = (int) a.size();
        std::vector<std::vector<int> > result(n, std::vector<int>(n, 0));
        long long int start = clock();
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < n; i++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        long long int end = clock();
        Log::print("kji", start, end);
    }
};

class Context {
protected:
    Strategy *operation;
public:
    virtual ~Context() {}

    virtual void execute(std::vector<std::vector<int> > a, std::vector<std::vector<int> > b) {
        operation->multiply(a, b);
    }

    virtual void setStrategy(Strategy *strategy) {
        operation = strategy;
    }
};

class Runner {
private:
    int size;
public:
    void setMatrixSize(int size = 0) {
        this->size = size;
    }

    void run() {
        MatrixGenerator matrixGenerator = MatrixGenerator(this->size);

        matrixGenerator.generate();
        std::vector<std::vector<int> > a = matrixGenerator.getGeneratedMatrix();

        matrixGenerator.generate();
        std::vector<std::vector<int> > b = matrixGenerator.getGeneratedMatrix();

        IJKStrategy ijkStrategy;

        IKJStrategy ikjStrategy;
        JIKStrategy jikStrategy;
        JKIStrategy jkiStrategy;
        KIJStrategy kijStrategy;
        KJIStrategy kjiStrategy;

        printf("Size of the matrix: %dx%d\n", this->size, this->size);

        Context context;

        context.setStrategy(&ijkStrategy);
        context.execute(a, b);

        context.setStrategy(&ikjStrategy);
        context.execute(a, b);

        context.setStrategy(&jikStrategy);
        context.execute(a, b);

        context.setStrategy(&jkiStrategy);
        context.execute(a, b);

        context.setStrategy(&kijStrategy);
        context.execute(a, b);

        context.setStrategy(&kjiStrategy);
        context.execute(a, b);

        printf("\n");
    }
};

int main(int argc, char const *argv[]) {
    Runner runner;
    for (int n = 500; n <= 5000; n += 500) {
        runner.setMatrixSize(n);
        runner.run();
    }

    for (int n = 6000; n <= 10000; n += 1000) {
        runner.setMatrixSize(n);
        runner.run();
    }

    return 0;
}