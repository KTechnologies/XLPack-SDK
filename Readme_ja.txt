XLPack SDK (XLPack 7.0�p) (2024.9.5)

1. XLPack SDK�Ƃ�

XLPack���g�p���Ă��鐔�l�v�Z�֐����W���[����, C/C++�Ȃ�Excel VBA�ȊO�̊J������
��g�p���邽�߂̃��C�u��������юg�p�Ⴉ��Ȃ�܂�. �{SDK���g�p���č쐬�����v��
�O������, XLPack���C���X�g�[�������p�\�R���Ŏ��s���邱�Ƃ��ł��܂�.

�� - �g�p���XLPack�̊�{�@�\(�A�h�I���Ȃ�)�͈̔͂ō쐬����Ă���, XLPack���C��
     �X�g�[������Ă���Ύ��s���邱�Ƃ��ł��܂�.

2. �{�\�t�g�E�F�A�̎g�p����

���g�p�̍ۂɂ͉��L�g�p�����ɏ]���Ă�������.

�E�{�\�t�g�E�F�A�Ɋւ���S�Ă̌����̓P�C�e�N�m���W�[�ɋA�����܂�.
�E�{�\�t�g�E�F�A����ъ֘A�Z�p����, XLPack�̓����C���^�[�t�F�[�X������̂܂܌�
  �J������̂�, ���ׂĂ̏������Ńe�X�g���ꂽ���̂ł͂���܂���. �]����, ������
  �@�\, ����, �M��������ۏ؂�����̂ł͂���܂���. �܂�, �\���Ȃ��d�l�ύX���s��
  ���Ƃ�����܂�.
�E�����p�ɂ������Ă͎g�p�҂̐ӔC�̂��Ƃł��g�p��������. �\�t�g�E�F�A�g�p�̌��ʐ�
  ���������Ȃ鑹�Q���ɂ��Ă��ۏ؂����܂���. �܂�, �{�\�t�g�E�F�A�ɂ̓T�|�[�g�T
  �[�r�X�͒񋟂���܂���.
�E�{�\�t�g�E�F�A���Ĕz�z���邱�Ƃ͂ł��܂���.

3. �{�\�t�g�E�F�A�̍\��

SDK�ɂ͎��̃t�@�C�����܂܂�Ă��܂�.

Readme_ja.txt - �����t�@�C�� (�{�t�@�C��)
Readme.txt - �����t�@�C�� (�p��)
include\
  cnumlib.h - C�C���^�[�t�F�[�X�E�w�b�_�[
  cnumlib_mangling.h - C�C���^�[�t�F�[�X�E�݊����w�b�_�[ (6���Q��)
  cnumlib_complex.h - C�C���^�[�t�F�[�X�E�w�b�_�[ (���f���^��`)
  cblas.h - CBLAS�w�b�_�[
  lapacke.h - LAPACKE�w�b�_�[
  spcnumlib.h - C�C���^�[�t�F�[�X�E�w�b�_�[ (�a�s��)
  spcnumlib_mangling.h - C�C���^�[�t�F�[�X�E�݊����w�b�_�[
  pdecnumlib.h - C�C���^�[�t�F�[�X�E�w�b�_�[ (�Δ���������)
  pdecnumlib_mangling.h - C�C���^�[�t�F�[�X�E�݊����w�b�_�[
lib\
  XLPack.lib - C�C���^�[�t�F�[�X�E���C�u���� (64�r�b�g��)
  XLPack_32.lib - C�C���^�[�t�F�[�X�E���C�u���� (32�r�b�g��)
  Lapacke.lib - LAPACKE/CBLAS���C�u���� (64�r�b�g��)
  Lapacke_32.lib - LAPACKE/CBLAS���C�u���� (32�r�b�g��)
samples\
  �g�p��: txt�t�@�C���̓R�}���h���C���Ŏ��s�����o�͗�ł�.
    C_C++ - C/C++�ɂ��g�p��
      testcpp_Matrix�́u�X�g���E�X�g���b�v�̃v���O���~���O����v2011, �ĉj�� ��
      �o�Ă���Matrix.h�����MatrixIO.h��ʓr�p�ӂ���K�v������܂�.
    Python - Python�ɂ��g�p��
      numpy���C���X�g�[������Ă���K�v������܂�.
      pyd�t�@�C����Windows�p��64�r�b�g�o�C�i���t�@�C���� Python 3.7 �ȍ~�Ŏg�p
      �ł��܂�. 32�r�b�gPython�̏ꍇ, ���邢�͂��܂����삵�Ȃ��ꍇ�ɂ�
      XLPack.py(�����̓��������ctypes��)���g�p���Ă�������.
    C# - C#�ɂ��g�p��

4. �����

���L�\�t�g�E�F�A���g�p���Ċm�F���s���Ă��܂�. �C���X�g�[����Ԃɂ���Ă͎g�p�Ⴛ
�̂܂܂ł͓��삹���C�����K�v�ȏꍇ������܂�.

�EXLPack 7.0.0
�EWindows 11 (23H2)
  �EVisual Studio 2022 (17.11.0)
  �EPython 3.12.5 (numpy 2.0.1)

5. �h�L�������g

�I�����C���}�j���A�������z�[���y�[�W(*)�Ō��J���Ă��܂�.

(*) https://www.ktech.biz/jp/

6. C�C���^�[�t�F�[�X�̕ύX

�{�o�[�W�����ł� XLPack 6.0 �܂ł� C/C++ �֐����̐擪�ɃA���_�[���C��(_)���ǉ���
��Ă��܂�. cnumlib_mangling.h �w�b�_�[�ɂ��݊����̂���֐�������`����܂�.

---
(C) 2014-2024  K Technologies
