
#ifndef  IO_utils_h
#define  IO_utils_h

// A simple function that will read a file into an allocated char pointer buffer
char* filetobuf(char const  *file){
	FILE *fptr;
	long length;
	char *buf;

	fptr = fopen(file, "rb");			// Open file for reading
	if (!fptr)							// Return NULL on failure
	    return NULL;
	fseek(fptr, 0, SEEK_END); 			// Seek to the end of the file
	length = ftell(fptr); 				// Find out how many bytes into the file we are
	buf = (char*)malloc(length+1); 		// Allocate a buffer for the entire length of the file and a null terminator
	fseek(fptr, 0, SEEK_SET); 			// Go back to the beginning of the file
	fread(buf, length, 1, fptr); 		// Read the contents of the file in to the buffer
	fclose(fptr); 						// Close the file
	buf[length] = 0; 					// Null terminator

	return buf; 						// Return the buffer
}

#endif
