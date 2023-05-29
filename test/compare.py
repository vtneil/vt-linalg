import subprocess

executables = ["./naive", "./strassen"]
n_values = list(range(1, 13))

print("N,naive,strassen")

for n in n_values:
    times = []
    for executable in executables:
        command = f"time ./{executable} {n}"
        try:
            t_acc = 0
            for i in range(5):
                output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
                output = output.decode("utf-8").strip()
                tim = output.splitlines()[0][5:].split("m")
                tim_n = int(tim[0]) * 60 + float(tim[1][:-1])
                t_acc += tim_n
            times.append(f"{t_acc / 10:.3f}")
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while executing {command}: {e.output.decode('utf-8').strip()}")
            times.append("N/A")
    print(f"{2 ** n},{times[0]},{times[1]}")
