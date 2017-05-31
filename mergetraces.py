import os

TIME_LIMIT = 30

def last_value_up_to(x, pairs):
    val = -1
    for pair in pairs:
        if pair[0] <= x:
            val = pair[1]
        else:
            break
    return val


def build_csv_from_results(merged_results):
    ostr = 't;' + ';'.join(trace_files) + '\n'
    for t in sorted(merged_results):
        ostr += str(t) + ';' + ';'.join(map(lambda n: str(n), merged_results[t])) + '\n'
    return ostr.replace('.',',')


trace_files = list(filter(lambda fn: fn.endswith('Trace.txt'), os.listdir('.')))

results_for_tfile = {}

for trace_file in trace_files:
    with open(trace_file) as fp:
        lines = fp.readlines()[1:]

        last_t = -1
        last_obj = -1

        results_for_tfile[trace_file] = []

        for line in lines:
            parts = line.split(';')
            t = float(parts[0])
            obj = float(parts[1])

            if t == last_t and obj > last_obj:
                results_for_tfile[trace_file][-1] = (t, obj)
                last_obj = obj
            elif t > last_t:
                results_for_tfile[trace_file].append((t, obj))
                last_obj = obj
                last_t = t

trange = [round(i * 0.01, 2) for i in range(100)] + [i + 1 for i in range(TIME_LIMIT)]
merged_results = {t: [last_value_up_to(t, results_for_tfile[trace_file]) for trace_file in trace_files] for t in trange}

with open('mergedtraces.txt', 'w') as fp:
    fp.write(build_csv_from_results(merged_results))
