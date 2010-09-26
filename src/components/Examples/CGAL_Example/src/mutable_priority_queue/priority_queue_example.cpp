#include<iostream>
#include"mutable_priority_queue.h"

template<typename Val>
class Compare {
public:
	const inline bool operator()( const Val a, const Val b ) const {
		return a < b;
	}
};

int main( int argc, char **argv ){

	Compare<int> c;
	mutable_priority_queue<float,int,std::greater<float>,Compare<int> > A(c), B;
	mutable_priority_queue<float,int,std::greater<float>,Compare<int> >::iterator i;

	A.insert( 1,  0.0f );
	A.insert( 2,  2.0f );
	A.insert( 3, -1.0f );

	std::cout << "A:" << std::endl;
	for( i=A.begin(); i!=A.end(); i++ ){
		std::cout << "key: " << i->first << ", value: " << i->second << std::endl;
	}

	A.update( 2, -5.0f );
	std::cout << "A:" << std::endl;
	for( i=A.begin(); i!=A.end(); i++ ){
		std::cout << "key: " << i->first << ", value: " << i->second << std::endl;
	}

	i = A.begin();

	// Can't do this (by design), gives a compile error
	//A.erase(i);

	// Do this instead
	A.erase( i->second );

	// after performing an erase, the iterators are 
	// invalid, so this gives a runtime error:
	// i++;

	std::cout << "A:" << std::endl;
	for( i=A.begin(); i!=A.end(); i++ ){
		std::cout << "key: " << i->first << ", value: " << i->second << std::endl;
	}

	std::cout << "B:" << std::endl;
	for( i=B.begin(); i!=B.end(); i++ ){
		std::cout << "key: " << i->first << ", value: " << i->second << std::endl;
	}
	B = A;
	B.update( 2, 5.0f );
	std::cout << "B:" << std::endl;
	for( i=B.begin(); i!=B.end(); i++ ){
		std::cout << "key: " << i->first << ", value: " << i->second << std::endl;
	}
	std::cout << "A:" << std::endl;
	for( i=A.begin(); i!=A.end(); i++ ){
		std::cout << "key: " << i->first << ", value: " << i->second << std::endl;
	}

	return 0;
};
