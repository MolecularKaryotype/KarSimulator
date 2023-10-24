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

x = "(yy)(zz)"
y = ""
print(y.split('('))
