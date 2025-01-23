﻿using CodeHelper;
using Controls;
using HardWares.温度控制器;
using HardWares.温度控制器.SRS_PTC10;
using HardWares.相机_CCD_;
using HardWares.端口基类;
using HardWares.端口基类部分;
using HardWares.纳米位移台;
using HardWares.纳米位移台.PI;
using ODMR_Lab.Windows;
using ODMR_Lab.基本窗口;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using 温度监控程序.Windows;
using ContextMenu = Controls.ContextMenu;

namespace ODMR_Lab.相机
{
    /// <summary>
    /// Page1.xaml 的交互逻辑
    /// </summary>
    public partial class DevicePage : PageBase
    {
        public List<CameraInfo> Cameras { get; set; } = new List<CameraInfo>();

        public DevicePage()
        {
            InitializeComponent();
        }

        public override void Init()
        {
        }

        public override void CloseBehaviour()
        {
        }

        private void NewConnect(object sender, RoutedEventArgs e)
        {
            ConnectWindow window = new ConnectWindow(typeof(CameraBase), MainWindow.Handle);
            bool res = window.ShowDialog(out PortObject dev);
            if (res == true)
            {
                CameraInfo camera = (CameraInfo)new CameraInfo().CreateDeviceInfo(dev as CameraBase, window.ConnectInfo);
                Cameras.Add(camera);
                RefreshPanels();
            }
            else
            {
                return;
            }
        }

        public void RefreshPanels()
        {
            DeviceList.ClearItems();
            foreach (var item in Cameras)
            {
                DeviceList.AddItem(item, item.Device.ProductName);
            }
        }

        /// <summary>
        /// 显示摄像头窗口
        /// </summary>
        /// <param name="arg1"></param>
        /// <param name="arg2"></param>
        private void ShowCameraWindow(int arg1, object arg2)
        {
            CameraInfo info = arg2 as CameraInfo;
            if (info.DisplayWindow != null)
            {
                info.DisplayWindow.Topmost = true;
                info.DisplayWindow.Topmost = false;
            }
        }

        /// <summary>
        /// 右键菜单事件
        /// </summary>
        /// <param name="arg1"></param>
        /// <param name="arg2"></param>
        /// <param name="arg3"></param>
        private void ContextMenuEvent(int arg1, int arg2, object arg3)
        {
            CameraInfo inf = arg3 as CameraInfo;
            #region 关闭设备
            if (arg1 == 0)
            {
                if (MessageWindow.ShowMessageBox("提示", "确定要关闭此设备吗？", MessageBoxButton.YesNo, owner: MainWindow.Handle) == MessageBoxResult.Yes)
                {
                    inf.CloseDeviceInfoAndSaveParams(out bool result);
                    if (result == false) return;
                    Cameras.Remove(inf);
                    RefreshPanels();
                }
            }
            #endregion

            #region 在独立窗口中显示
            if (arg1 == 1)
            {
                if (inf.DisplayWindow != null)
                {
                    inf.DisplayWindow.Show();
                }
                else
                {
                    CameraWindow window = new CameraWindow(inf);
                    inf.DisplayWindow = window;
                    window.Show();
                }
            }
            #endregion

            #region 参数设置
            if (arg1 == 2)
            {
                ParameterWindow window = new ParameterWindow(inf.Device.AvailableParameterNames());
                window.ShowDialog();
            }
            #endregion
        }
    }
}
