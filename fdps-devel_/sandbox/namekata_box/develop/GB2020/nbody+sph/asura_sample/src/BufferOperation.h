/* BufferOperation.c */
void InitializeCommunicationBuffers(void);
void ReportAllocatedMemorySizes(void);
void CheckSizeofBufferExportSend(const int NumberofExportSend[restrict], const size_t SizeofElement);
void CheckSizeofBufferExportSendIndex(const int NumberofExportSend, const size_t SizeofElement, const int Index);
void CheckSizeofBufferExportRecv(const int NumberofExportRecv, const size_t SizeofElement);
void CheckSizeofBufferImportSend(const int NumberofImportSend, const size_t SizeofElement);
void CheckSizeofBufferImportRecv(const int NumberofImportRecv[restrict], const size_t SizeofElement);
void CheckSizeofBufferImportRecvIndex(const int NumberofImportRecv, const size_t SizeofElement, const int Index);
void ReleaseAllBuffers(void);
