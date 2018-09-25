///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAMEMORY_HPP
#define INAMEMORY_HPP

#include <cstdint>

/**
 * Proposes methods to allocate aligned memory.
 * Any memory from this class must be dellocated with
 * the Delete method.
 * Note: the alignement must be a power of 2
 */
class InaMemory {
    const std::size_t HeaderSize = sizeof(unsigned char*) + sizeof(std::size_t);

    template < std::size_t Alignement >
    static void* Allocate(const std::size_t inSize, const std::size_t inHeader) {
        // Return null for empty allocation
        if (inSize == 0) {
            return nullptr;
        }

        // Ensure it is a power of 2 > 0
        static_assert(Alignement != 0 && ((Alignement - 1) & Alignement) == 0, "Alignement must be a power of 2");

        // We will need to store the adress of the real blocks
        const std::size_t Offset = (Alignement < HeaderSize ? HeaderSize / Alignement : 1) * Alignement;

        unsigned char* allocatedMemoryPtr = new unsigned char[inSize + Offset - 1];
        unsigned char* alignedMemoryPtr   = reinterpret_cast< unsigned char* >((reinterpret_cast< std::size_t >(allocatedMemoryPtr) + Offset - 1) & ~(Alignement - 1));
        unsigned char* headerPtr          = (alignedMemoryPtr - HeaderSize);

        // Save real address then header value then alignement
        *reinterpret_cast< unsigned char** >(headerPtr)                       = allocatedMemory;
        *reinterpret_cast< std::size_t* >(headerPtr + sizeof(unsigned char*)) = inHeader;

        // Return aligned address
        return reinterpret_cast< void* >(alignedMemoryAddress);
    }

    static void Release(const void* ptrToFree) {
        // If equal to null do nothing
        if (ptrToFree) {
            // Retrieve real address
            const unsigned char* const* storeRealAddress = reinterpret_cast< const unsigned char* const* >(reinterpret_cast< const unsigned char* >(ptrToFree) - HeaderSize);
            delete[] reinterpret_cast< const unsigned char* >(*storeRealAddress);
        }
    }

public:
    template < std::size_t Alignement, class ObjectType >
    static ObjectType* _new(const std::size_t inNbElementsInArray) {
        if (inNbElementsInArray == 0) {
            return nullptr;
        }

        const std::size_t sizeInBytes = (inNbElementsInArray * sizeof(ObjectType));

        ObjectType* alignedArray = reinterpret_cast< ObjectType* >(Allocate(sizeInBytes, inNbElementsInArray));

        for (std::size_t idx = 0; idx < inNbElementsInArray; ++idx) {
            new (&alignedArray[idx]) ObjectType();
        }
        return array;
    }

    template < class ObjectType >
    inline void _delete(const ObjectType* ptrToFree) {
        if (ptrToFree) {
            const std::size_t numberOfElements = (*reinterpret_cast< const std::size_t* >(reinterpret_cast< const unsigned char* >(ptrToFree) + sizeof(unsigned char*) - HeaderSize));

            for (std::size_t idx = 0; idx < numberOfElements; ++idx) {
                ptrToFree[idx].~ObjectType();
            }

            Release(ptrToFree);
        }
    }
};


#endif // INAMEMORY_HPP
