﻿using CodeHelper;
using Controls;
using Controls.Windows;
using HardWares.Lock_In;
using HardWares.Windows;
using HardWares.仪器列表.电动翻转座;
using HardWares.射频源;
using HardWares.温度控制器;
using HardWares.温度控制器.SRS_PTC10;
using HardWares.源表;
using HardWares.相机_CCD_;
using HardWares.端口基类;
using HardWares.端口基类部分;
using HardWares.纳米位移台;
using HardWares.纳米位移台.PI;
using ODMR_Lab.Windows;
using ODMR_Lab.基本窗口;
using ODMR_Lab.设备部分;
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
using ContextMenu = Controls.ContextMenu;

namespace ODMR_Lab.设备部分.射频源_锁相放大器
{
    /// <summary>
    /// Page1.xaml 的交互逻辑
    /// </summary>
    public partial class DevicePage : DevicePageBase
    {
        public override string PageName { get; set; } = "信号发生器/锁相放大器";
        public List<SignalGeneratorInfo> SignalGenerators { get; set; } = new List<SignalGeneratorInfo>();
        public List<LockinInfo> LockIns { get; set; } = new List<LockinInfo>();


        public DevicePage()
        {
            InitializeComponent();
        }

        public override void InnerInit()
        {
        }

        public override void CloseBehaviour()
        {
        }

        public override void UpdateParam()
        {
        }

        private void NewRFSourceConnect(object sender, RoutedEventArgs e)
        {
            ConnectWindow window = new ConnectWindow(typeof(SignalGeneratorBase));
            bool res = window.ShowDialog(Window.GetWindow(this));
            if (res == true)
            {
                SignalGeneratorInfo rfsource = new SignalGeneratorInfo() { Device = window.ConnectedDevice as SignalGeneratorBase, ConnectInfo = window.ConnectInfo };
                rfsource.CreateDeviceInfoBehaviour();

                SignalGenerators.Add(rfsource);
                RefreshPanels();
            }
            else
            {
                return;
            }
        }
        private void NewLockInConnect(object sender, RoutedEventArgs e)
        {
            ConnectWindow window = new ConnectWindow(typeof(LockInBase));
            bool res = window.ShowDialog(Window.GetWindow(this));
            if (res == true)
            {
                LockinInfo lockin = new LockinInfo() { Device = window.ConnectedDevice as LockInBase, ConnectInfo = window.ConnectInfo };
                lockin.CreateDeviceInfoBehaviour();

                LockIns.Add(lockin);
                RefreshPanels();
            }
            else
            {
                return;
            }
        }

        public override void RefreshPanels()
        {
            SignalGeneratorList.ClearItems();
            LockInList.ClearItems();

            foreach (var item in SignalGenerators)
            {
                SignalGeneratorList.AddItem(item, item.Device.ProductName);
            }

            foreach (var item in LockIns)
            {
                LockInList.AddItem(item, item.Device.ProductName);
            }
        }

        /// <summary>
        /// 右键菜单事件
        /// </summary>
        /// <param name="arg1"></param>
        /// <param name="arg2"></param>
        /// <param name="arg3"></param>
        private void SignalGeneratorContextMenuEvent(int arg1, int arg2, object arg3)
        {
            SignalGeneratorInfo inf = arg3 as SignalGeneratorInfo;
            #region 关闭设备
            if (arg1 == 0)
            {
                if (MessageWindow.ShowMessageBox("提示", "确定要关闭此设备吗？", MessageBoxButton.YesNo, owner: Window.GetWindow(this)) == MessageBoxResult.Yes)
                {
                    inf.CloseDeviceInfoAndSaveParams(out bool result);
                    if (result == false) return;
                    SignalGenerators.Remove(inf);
                    ChannelList.ClearItems();
                    RefreshPanels();
                }
            }
            #endregion

            #region 参数设置
            if (arg1 == 1)
            {
                ParameterWindow window = new ParameterWindow(inf.Device, Window.GetWindow(this));
                window.ShowDialog();
            }
            #endregion
        }


        private void LockInContextMenuEvent(int arg1, int arg2, object arg3)
        {
            LockinInfo inf = arg3 as LockinInfo;
            #region 关闭设备
            if (arg1 == 0)
            {
                if (MessageWindow.ShowMessageBox("提示", "确定要关闭此设备吗？", MessageBoxButton.YesNo, owner: Window.GetWindow(this)) == MessageBoxResult.Yes)
                {
                    inf.CloseDeviceInfoAndSaveParams(out bool result);
                    if (result == false) return;
                    LockIns.Remove(inf);
                    RefreshPanels();
                }
            }
            #endregion

            #region 参数设置
            if (arg1 == 1)
            {
                ParameterWindow window = new ParameterWindow(inf.Device, Window.GetWindow(this));
                window.ShowDialog();
            }
            #endregion
        }

        private void SignalGeneratorSelectEvent(int arg1, object arg2)
        {
            //选中对应发生器时刷新通道列表
            ChannelList.ClearItems();
            SignalGeneratorInfo inf = arg2 as SignalGeneratorInfo;
            foreach (var item in inf.Channels)
            {
                ChannelList.AddItem(item, inf.Device.ProductName + " " + item.Device.ChannelName, item.Device.IsOutOpen);
            }
        }

        private void ChannelContextMenuEvent(int arg1, int arg2, object arg3)
        {
            SignalGeneratorChannelInfo inf = arg3 as SignalGeneratorChannelInfo;
            #region 参数设置
            if (arg1 == 0)
            {
                ParameterWindow window = new ParameterWindow(inf.Device, Window.GetWindow(this));
                window.ShowDialog();
            }
            #endregion
        }
    }
}
