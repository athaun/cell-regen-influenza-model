# from datetime import datetime
# dateNow = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
# print('\n[{}] INFO: A cell regeneration has occured\n'.format(dateNow))

def progbar(curr, total, full_progbar):
    frac = curr/total
    filled_progbar = round(frac*full_progbar)
    print('\r', '#'*filled_progbar + '-'*(full_progbar-filled_progbar), '[{:>7.2%}]'.format(frac), end='')

for i in range(100000+1):
    progbar(i, 100000, 20)
print()
