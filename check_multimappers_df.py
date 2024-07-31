reads_count = 100
user_batch_size = 0
default_batch_size = 50

# calculate batch size
if (user_batch_size and (user_batch_size>0)):
    batch_size = user_batch_size
elif (reads_count and (reads_count>0)):
    batch_size = int(reads_count/1)
else:
    batch_size =default_batch_size

print(batch_size)