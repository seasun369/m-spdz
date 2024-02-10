def format_and_write_data(input_filename, output_filename):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        for line in infile:
            # 移除行尾的换行符，然后以空格分割行内数据
            parts = line.strip().split()
            row_number = parts[0]  # 提取行号
            data_parts = parts[1:]  # 提取除行号外的数据部分
            # 将数据之间的空格替换为半角逗号
            formatted_data = ', '.join(data_parts)
            # 在行号之后立即加上左大括号，并在整行数据的末尾加上右大括号
            formatted_line = f"{row_number} {{{formatted_data}}}"
            # 写入到输出文件，末尾添加换行符
            outfile.write(formatted_line + '\n')


input_filename = 'mac_c.txt'
output_filename = 'direct_mac_c.txt'
format_and_write_data(input_filename, output_filename)
