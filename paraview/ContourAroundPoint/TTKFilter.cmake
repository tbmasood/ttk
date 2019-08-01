# Allows to disable each filter
option(TTK_BUILD_CONTOUR_AROUND_POINT_FILTER "Build the ContourAroundPoint filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_CONTOUR_AROUND_POINT_FILTER)

if(${TTK_BUILD_CONTOUR_AROUND_POINT_FILTER})
  ttk_register_pv_filter(ttkContourAroundPoint ContourAroundPoint.xml)
endif()

