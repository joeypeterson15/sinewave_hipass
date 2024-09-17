#include <iostream>
#include "kiss_fftr.h"
#include <cmath>
#include <fstream>
#include <array>
#include <string>
#include <complex>
#include <vector>
#include <iomanip>
#include <cstdint>

using namespace std;

const int sampleRate = 44100;
const int bitDepth = 16;
const float nyquisFrequency = sampleRate / 2;
const float c4NoteHz = 262;
const float d4Hz = 293.66;
const float e5NoteHz = 659.25;
const float f4Hz = 349.23;
const float g5NoteHz = 783.99;
const float a4NoteHz = 440;
const float b4Hz = 493.88;
const double e = 2.71828182845904523536;


class SineOscillator {
    float frequency, amplitude, angle = 0.0f, offset = 0.0f;
public:
    SineOscillator(float freq, float amp) : frequency(freq), amplitude(amp) {
        offset = 2 * M_PI * frequency / sampleRate;
    }
    float process() {
        auto sample = amplitude * sin(angle);
        angle += offset;
        return sample; 
        // Asin(2(pi)f/samplerate)
    }
};

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

void writeToFile(ofstream &file, int value, int size) {
    file.write(reinterpret_cast<const char*> (&value), size);
}

void generateWaveFile(int N, std::vector<float> &data, string fileSuffix) {
    ofstream audioFile;
    string fileName = "Amin7" + fileSuffix + ".wav";
    audioFile.open(fileName, ios::binary);
    //Header chunk
    audioFile << "RIFF";
    audioFile << "----";
    audioFile << "WAVE";

    // Format chunk
    audioFile << "fmt ";
    writeToFile(audioFile, 16, 4); // Size
    writeToFile(audioFile, 1, 2); // Compression code
    writeToFile(audioFile, 1, 2); // Number of channels
    writeToFile(audioFile, sampleRate, 4); // Sample rate
    writeToFile(audioFile, sampleRate * bitDepth / 8, 4 ); // Byte rate
    writeToFile(audioFile, bitDepth / 8, 2); // Block align
    writeToFile(audioFile, bitDepth, 2); // Bit depth

    //Data chunk
    audioFile << "data";
    audioFile << "----";

    int preAudioPosition = audioFile.tellp();

    auto maxAmplitude = pow(2, bitDepth - 1) - 1;
    int counter = 0;
    for(int i = 0; i < N; i++ ) {
        auto sample = data[i];
        int intSample = static_cast<int> (sample * maxAmplitude);
        writeToFile(audioFile, intSample, 2);
    }

    int postAudioPosition = audioFile.tellp();

    audioFile.seekp(preAudioPosition - 4);
    writeToFile(audioFile, postAudioPosition - preAudioPosition, 4);

    audioFile.seekp(4, ios::beg);
    writeToFile(audioFile, postAudioPosition - 8, 4);

    audioFile.close();
}

int waveFileReader() {
    std::ifstream file("/Users/joeypeterson/Desktop/test_impulse_response.wav", std::ios::binary);

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

int main(int argc, char **argv) {
    // read impulse response test file
    waveFileReader();

    // string chord = argv[1]; //7th chord desired
    int duration = std::stoi(argv[1]); //audio length(seconds)
    int hiPassFrequencyCutoff = std::stoi(argv[2]);
    int N = sampleRate * duration; //number of samples

    SineOscillator sineOscillatorA(a4NoteHz, 0.125);
    SineOscillator sineOscillatorC(c4NoteHz, 0.125);
    SineOscillator sineOscillatorE(e5NoteHz, 0.125);
    SineOscillator sineOscillatorG(g5NoteHz, 0.125);

    // collect sound data
    std::vector <float> samples(N, 0.0f);
    for (int i = 0; i < N; ++i) {
        float sample = sineOscillatorA.process() + sineOscillatorC.process() + sineOscillatorE.process() + sineOscillatorG.process();
        samples[i] = sample;
    }

    // perform FFT with kissfft library 
    kiss_fftr_cfg forward_cfg = kiss_fftr_alloc(N, 0, NULL, NULL);
    std::vector<kiss_fft_cpx> outFreq(N / 2 + 1);
    std::fill(outFreq.begin(), outFreq.end(), kiss_fft_cpx{0.0, 0.0}); //watching for any plan/cache

    kiss_fftr(forward_cfg, samples.data(), outFreq.data());

    // process HIGH pass filter
    std::vector<kiss_fft_cpx> hpSamples(N / 2 + 1);
    int frequencyCutOffBinK = hiPassFrequencyCutoff * sampleRate / (N / 2 + 1);
    std::cout << "frequency cutoff bin: " << frequencyCutOffBinK << std::endl;
    for (int k = frequencyCutOffBinK; k >= 0; --k) {
        outFreq[k].r = 0.0;
        outFreq[k].i = 0.0;
    }

    // convert data back to time domain
    kiss_fftr_cfg inverse_cfg = kiss_fftr_alloc(N / 2, 1, NULL, NULL);
    std::vector<float> outTime(N);
    std::fill(outTime.begin(), outTime.end(), float(0.0)); //watching for any plan/cache
    kiss_fftri(inverse_cfg, outFreq.data(), outTime.data());

    free(inverse_cfg);
    free(forward_cfg); // Free the Kiss FFT configuration

    for (int i = 0; i < N; i++) {
        outTime[i] = outTime[i] / float(N); //normalize data
    }

    generateWaveFile(N, outTime, argv[2]);

    return 0;
}