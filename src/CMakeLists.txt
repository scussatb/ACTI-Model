file(GLOB SRC
	"*.h"
	"*.hpp"
	"*.cpp"
	)
add_executable(mecaBase ${SRC})
target_link_libraries(mecaBase mecacell)
#target_link_libraries(mecaBase mecacellviewer)
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
#find_package(Qt5Core)
#find_package(Qt5Widgets)
#find_package(Qt5Gui)
#find_package(Qt5Quick)
#find_package(Qt5OpenGL)
#qt5_use_modules(mecaBase Quick Core Gui OpenGL)
