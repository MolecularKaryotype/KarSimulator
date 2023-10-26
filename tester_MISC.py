from Preparation.Get_Chr_Names import *

# def outer_f(param1, param2, param3):
#     print(param1, param2, param3)
#
#
# def inner_f(param1, param2):
#     return param1, param2
#
#
# outer_f(1, *inner_f(2, 3))

# Get_All_Chr_Names('./Precompiled_Kar/Down_Extra21q.fasta', 'fullset_names.txt')

# a = [1]
# b = [2, 3]
# c = [4, 5]
# print([*a, *b, *c])


# class MyException(Exception):
#     pass
#
#
# def f1(input_int):
#     if input_int == 1:
#         raise MyException
#     else:
#         print(input_int)
#
#
# x = [j for j in range(0, 4)]
# pointer = 0
# for i in range(0, 5):
#     executed = False
#     while not executed:
#         try:
#             f1(x[pointer])
#         except MyException:
#             pointer = (pointer + 1) % 4
#             continue
#         executed = True
#         pointer = (pointer + 1) % 4

# my_dict = {
#     "banana": 2,
#     "apple": 5,
#     "date": 1,
#     "cherry": 8
# }
# sorted_items = sorted(my_dict.items(), key=lambda x: x[1], reverse=True)
# print(sorted_items[0][0])

# x = "(yy)(zz)"
# y = ""
# print(y.split('('))

# # Define the source list (the list you want to insert)
# source_list = [7, 8, 9]
# # Define the destination list (the list into which you want to insert)
# destination_list = [1, 2, 3, 4, 5, 6]
# # Specify the index where you want to insert the source list
# insert_index = 3
# # Use slicing and the extend method to insert the source list into the destination list
# destination_list[insert_index:insert_index] = source_list
# # Now, destination_list will contain the elements from the source_list inserted at the specified index
# print(destination_list)

# left = 3
# right = 5
# for i in range(right-1,left-1,-1):
#     print(i)

x = '\t1\t2\t'
print(len(x.split('\t')))
print(x.split('\t'))
