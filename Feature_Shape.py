Featuresval = open("methylation_recovery_colnames_values_re.txt")
Featuresval = Featuresval.read()
Featuresval = Featuresval.strip().split("\n")
for i in Featuresval:
    i = i.split("\t")
    #print(i[1])
    #print(type(int(i[1])))
    FShaped = open("methylation_recovery_colnames_values_re_shaped.txt", "a")
    if int(i[1]) >= 0 and int(i[1]) <= 100:
        i.insert(2, "1\n")
        range0to100 = "\t".join(i)
        #print(range0to100)
        FShaped.write(range0to100)
        FShaped.close()
    elif int(i[1]) > 100 and int(i[1]) <= 500:
        i.insert(2, "2\n")
        range100to500 = "\t".join(i)
        #print(range100to500)
        FShaped.write(range100to500)
        FShaped.close()
    elif int(i[1]) > 500 and int(i[1]) <= 1000:
        i.insert(2, "3\n")
        range500to1000 = "\t".join(i)
        FShaped.write(range500to1000)
        FShaped.close()
    elif int(i[1]) > 1000 and int(i[1]) <= 10000:
        i.insert(2, "4\n")
        range1000to10000 = "\t".join(i)
        FShaped.write(range1000to10000)
        FShaped.close()
    elif int(i[1]) > 10000:
        i.insert(2, "5\n")
        rangeabv10000 = "\t".join(i)
        FShaped.write(rangeabv10000)
        FShaped.close()
    else:
        print("Negative value encountered")

print("Completed")
