{
  'targets': [
    {
      'target_name': 'mist_math',
      'type': 'shared_library',
      'dependencies': [
      ],
      'include_dirs': [
        '..',
      ],
      'sources': [
        'degree.cc',
        'degree.h',
        'math_util.cc',
        'math_util.h',
        'radian.cc',
        'radian.h',
        'vector2.cc',
        'vector2.h',
        'vector3.cc',
        'vector3.h',
        'vector4.cc',
        'vector4.h',
      ],
    },
  ], # targets
}
