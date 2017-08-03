io.write("Calling howdy() ...\n")
local value = howdy("First", "Second", 112233)
io.write(string.format("howdy() returned: %s\n", tostring(value)))
