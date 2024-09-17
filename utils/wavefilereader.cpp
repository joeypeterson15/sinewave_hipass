#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>
#include <iomanip>

#pragma pack(push, 1)  // Ensures the structure is tightly packed

struct WAVHeader {
    char chunkID[4];        // "RIFF"
    uint32_t chunkSize;     // Size of the file minus 8 bytes
    char format[4];         // "WAVE"
    char subchunk1ID[4];    // "fmt "
    uint32_t subchunk1Size; // Size of the fmt chunk
    uint16_t audioFormat;   // Audio format (1 for PCM)
    uint16_t numChannels;   // Number of channels
    uint32_t sampleRate;    // Sample rate
    uint32_t byteRate;      // Byte rate
    uint16_t blockAlign;    // Block align
    uint16_t bitsPerSample; // Bits per sample
    char subchunk2ID[4];    // "data"
    uint32_t subchunk2Size; // Size of the data section
};

#pragma pack(pop)  // Restores default packing

int main(std::string waveFile) {
    std::ifstream file(waveFile, std::ios::binary);

    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    WAVHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(header));  // Read the header

    // Output header information in a readable format
    std::cout << "ChunkID: " << std::string(header.chunkID, 4) << std::endl;
    std::cout << "ChunkSize: " << header.chunkSize << std::endl;
    std::cout << "Format: " << std::string(header.format, 4) << std::endl;
    std::cout << "Subchunk1ID: " << std::string(header.subchunk1ID, 4) << std::endl;
    std::cout << "AudioFormat: " << header.audioFormat << std::endl;
    std::cout << "NumChannels: " << header.numChannels << std::endl;
    std::cout << "SampleRate: " << header.sampleRate << std::endl;
    std::cout << "BitsPerSample: " << header.bitsPerSample << std::endl;
    std::cout << "Subchunk2Size: " << header.subchunk2Size << std::endl;

    // Optionally, read audio data (not human-readable) into a buffer
    std::vector<char> audioData(header.subchunk2Size);
    file.read(audioData.data(), header.subchunk2Size);

    file.close();
    return 0;
}
