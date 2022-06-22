from glob import glob
import os
import csv

id2problem = {
    '1': 'mktrap',
    '2': 'ftrap',
    '3': 'cyctrap',
    '4': 'nk',
    '5': 'spin',
    '6': 'sat',
    '7': 'maxcut',
    '8': 'USal_NSize_3_7',
    '9': 'USal_NSize_3_10',
    '10': 'linear_mktrap',
    '11': 'exponential_mktrap',
    '12': 'four_five_six',
    '13': 'three_four_five_six_seven',
    '14': 'ftrap4',
    '15': 'ftrap6'
}

rows = []

for dir in glob('./out_dir/*'):
    csv_list = glob(f'{dir}/*.csv')
    s_list = [x.split('/')[-1] for x in csv_list if len(x.split('/')[-1].split('_'))==2]
    for s in s_list:
        problem_str = s.split('_')
        problem_name = id2problem[problem_str[0]]
        problem_len = problem_str[1].split('.')[0]
        with open(os.path.join(dir, s), 'r') as f:
            lines = f.readlines()
        row1 = lines[1].split()
        row2 = lines[2].split()
        # print(dir.split('/')[-1], end=' ')
        # print(tmp[1], tmp[2])
        print([dir.split('/')[-1], problem_name, problem_len, row1[1], row2[1]])
        # rows.append([dir.split('/')[-1], problem_name, problem_len, row1[1], row2[1]])
        rows.append([dir.split('/')[-1], problem_name, problem_len, row1[0], row2[0], row1[1], row2[1], row1[4], row2[4], row1[5], row2[5]])


with open('summary.csv', 'w') as f:
      
    # using csv.writer method from CSV package
    write = csv.writer(f)
      
    write.writerow(['GA', 'Problem', 'ell', 'NFE_mean', 'NFE_std', 'RM_mean', 'RM_std', 'BM_mean', 'BM_std'])
    write.writerows(rows)
